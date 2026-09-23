import os,sys
import gzip
import argparse
import random
import pandas as pd
from reportlab.pdfgen import canvas
import PyPDF2
import pdfplumber
import re


# Check if the specified directory exists; create it if it doesn't
def create_directory_if_not_exists(directory):
    if not os.path.exists(directory):
        os.makedirs(directory, exist_ok=True)


# Extract the width and height of the first page of a PDF
def get_pdf_dimensions(pdf_file):
    reader = PyPDF2.PdfReader(pdf_file)
    first_page = reader.pages[0]
    pdf_width = float(first_page.mediabox.width)
    pdf_height = float(first_page.mediabox.height)
    return pdf_width, pdf_height


# Parse a GTF file and extract exons for each transcript within a specified genomic region.
# The function returns a dictionary with transcripts as keys and lists of exon details as values.
def read_gtf_and_number_exons(gtf_file, region):
    """
    Reads a GTF file and extracts exon information grouped by transcript.

    Args:
        gtf_file (str): Path to the GTF file.
        region (str): Genomic region in "chr:start-end" format.

    Returns:
        dict: Dictionary {transcript_id: [{"start": int, "end": int, "exon_number": int}]}
    """
    
    # Extract chromosome and coordinates from the region string
    chr_region = region.split(":")[0]
    start = int(region.split(":")[1].split("-")[0])
    end = int(region.split(":")[1].split("-")[1])
    exons = {}

    # Open the GTF file (handle gzip compression if needed)
    with gzip.open(gtf_file, 'rt') if gtf_file.endswith(".gz") else open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):  # Skip comment lines
                continue
            
            # Split the GTF fields
            fields = line.strip().split("\t")
            gtf_chr, _, feature_type, el_start, el_end, _, strand, _, attributes = fields

            # Ignore features that are not exons or outside the region
            if gtf_chr != chr_region or feature_type != "exon":
                continue

            # Convert coordinates to 0-based indexing
            el_start, el_end = int(el_start) - 1, int(el_end)

            # Extract transcript ID using regex
            try:
                transcript_id = re.search(r'transcript_id "([^"]+)"', attributes).group(1)
            except AttributeError:
                print(f"Error: 'transcript_id' not found in line: {line}")
                continue

            # Extract exon number, if available
            exon_number = None
            if "exon_number" in attributes:
                try:
                    exon_number = int(re.search(r'exon_number (\d+)', attributes).group(1))
                except AttributeError:
                    print(f"Warning: 'exon_number' not found in line: {line}")

            # Store exons that overlap with the region of interest
            if el_end > start and el_start < end:
                exon_data = {
                    "start": max(el_start, start),
                    "end": min(el_end, end),
                    "exon_number": exon_number
                }
                exons.setdefault(transcript_id, []).append(exon_data)

    # Sort exons by genomic start position for each transcript
    for transcript_id in exons:
        exons[transcript_id].sort(key=lambda x: x["start"])

    return exons


# Check if a given Y position in the PDF is already occupied by text
def is_position_occupied(y, words, pdf_height):
    for word in words:
        word_top = pdf_height - word["top"]
        word_bottom = pdf_height - word["bottom"]
        if y >= word_bottom and y <= word_top:
            return True
    return False


# Compute transcript positions in a PDF and generate exon sequence labels.
# Adjusts label positioning based on available space and prevents overlap.
def calculate_transcript_labels(pdf_file, exon_annotations, pdf_width, pdf_height, nm, margin=1, gtf_file=None):
    """
    Calculates Y positions for each transcript and generates a label 
    containing the sequence of exon numbers.

    Args:
        pdf_file (str): Path to the PDF file (for extracting Y coordinates).
        exon_annotations (dict): Dictionary containing exon details for each transcript.
        pdf_width (float): PDF width.
        pdf_height (float): PDF height.
        margin (float): Left margin for positioning labels.
        nm (str): Transcript of reference (for coloring).
        gtf_file (str, optional): Path to GTF file for fallback.

    Returns:
        dict: Dictionary {transcript_id: {"x": float, "y": float, "label": str, "color": bool}}.
    """
    positions_by_transcript = {}

    with pdfplumber.open(pdf_file) as pdf:
        page = pdf.pages[0]
        words = page.extract_words()

        # Parcourir tous les transcripts visibles dans le PDF
        for word in words:
            transcript = word["text"]
            if not (transcript.startswith("NM_") or transcript.startswith("ENST")):
                continue  # ignorer les autres mots

            transcript_x = word["x0"]
            transcript_y = pdf_height - word["top"]

            # Cas 1 : exons connus
            exons = exon_annotations.get(transcript, [])
            if exons:
                exon_numbers = [str(exon["exon_number"]) for exon in exons]
                label = "-".join(exon_numbers)

                # === Ajustement largeur du label ===
                canvas_temp = canvas.Canvas(None)
                canvas_temp.setFont("Helvetica", 10)
                x_start = margin
                x_end = transcript_x
                max_label_width = abs(x_end - x_start)

                while (canvas_temp.stringWidth(label, "Helvetica", 10) > max_label_width 
                       and len(exon_numbers) > 2):
                    middle_index = len(exon_numbers) // 2

                    if len(exon_numbers) % 2 != 0:  # impair
                        exon_numbers.pop(middle_index)
                    else:  # pair
                        left_exon = int(exon_numbers[middle_index - 1])
                        right_exon = int(exon_numbers[middle_index])
                        if right_exon > left_exon:
                            exon_numbers.pop(middle_index)
                        else:
                            exon_numbers.pop(middle_index - 1)

                    if len(exon_numbers) == 2:
                        label = "".join(exon_numbers[:middle_index] + ["..."] + exon_numbers[middle_index:])
                    else:
                        label = "-".join(exon_numbers[:middle_index] + ["..."] + exon_numbers[middle_index:])

            # Cas 2/3 : aucun exon → fallback GTF
            else:
                label = None
                if gtf_file:
                    leftmost = None
                    with gzip.open(gtf_file, 'rt') if gtf_file.endswith(".gz") else open(gtf_file) as f:
                        for line in f:
                            if line.startswith("#"):
                                continue
                            fields = line.strip().split("\t")
                            if len(fields) < 9:
                                continue
                            gtf_chr, _, feature_type, el_start, el_end, _, _, _, attributes = fields
                            if feature_type != "exon":
                                continue
                            if f'transcript_id "{transcript}' not in attributes:
                                continue
                            el_start, el_end = int(el_start) - 1, int(el_end)
                            exon_number = None
                            if "exon_number" in attributes:
                                m = re.search(r'exon_number (\d+)', attributes)
                                if m:
                                    exon_number = int(m.group(1))
                            if leftmost is None or el_start < leftmost["start"]:
                                leftmost = {"start": el_start, "end": el_end, "exon_number": exon_number}
                    if leftmost and leftmost["exon_number"]:
                        label = str(leftmost["exon_number"])

            # Stockage
            positions_by_transcript[transcript] = {
                "x": transcript_x,
                "y": transcript_y - 10,
                "label": label,
                "color": nm in transcript
            }

        # Calcul du max_y pour éviter chevauchements
        if positions_by_transcript:
            max_y = max(data["y"] for data in positions_by_transcript.values())
            y_increment = 1
            while is_position_occupied(max_y, words, pdf_height):
                max_y += y_increment
        else:
            max_y = 0

    return positions_by_transcript, max_y


# Generate a PDF overlay containing transcript labels and genomic information.
# The labels are aligned to the left and positioned dynamically based on available space.
def create_labels_overlay(output_overlay, positions, pdf_width, pdf_height, chr, start, end, gene, event_type, max_y, margin=1):
    """
    Creates a PDF overlay with transcript labels aligned to the left.

    Args:
        output_overlay (str): Path to the generated overlay file.
        positions (dict): Dictionary {transcript_id: {"x": float, "y": float, "label": str}}.
        pdf_width (float): PDF width.
        pdf_height (float): PDF height.
        chr (str): Chromosome identifier.
        start (int): Start position of the genomic region.
        end (int): End position of the genomic region.
        gene (str): Gene associated with the event.
        event_type (str): Type of genomic event.
        max_y (float): Maximum Y position to align labels.
        margin (float): Left margin for positioning labels.

    Returns:
        None: Saves the generated overlay PDF.
    """
    
    # Initialize PDF canvas
    c = canvas.Canvas(output_overlay, pagesize=(pdf_width, pdf_height))
    c.setFont("Helvetica", 10)

    # Convert genomic information into formatted strings
    chr = f"chr : {chr}"
    gene = f"gene : {gene}"
    event_type = f"event : {event_type}"
    start = f"start : {start}"
    stop = f"end : {end}"

    # Draw genomic details at the top of the PDF overlay
    c.drawString(margin, max_y + 125, chr)
    c.drawString(margin, max_y + 105, gene)
    c.drawString(margin, max_y + 85, event_type)
    c.drawString(margin, max_y + 65, start)
    c.drawString(margin, max_y + 45, stop)

    # Label for exon numbering
    c.drawString(margin, max_y + 25, "exon number")

    # Iterate over each transcript and draw its label on the overlay
    for transcript, data in positions.items():
        label = data["label"]
        x_pos = margin  # Left alignment for labels
        y_pos = data["y"]
        color = data["color"]
        if color:
            c.setFillColorRGB(1, 0, 0)
        else:
            c.setFillColorRGB(0, 0, 0)
        c.drawString(x_pos, y_pos, label)

    # Save the PDF overlay
    c.save()


# Merge an annotation overlay PDF with the base Sashimi plot PDF.
# This ensures that transcript labels and genomic information are incorporated into the final output.
def merge_pdfs(base_pdf, overlay_pdf, output_pdf):
    with open(base_pdf, "rb") as base_file, open(overlay_pdf, "rb") as overlay_file:
        base_reader = PyPDF2.PdfReader(base_file)
        overlay_reader = PyPDF2.PdfReader(overlay_file)
        writer = PyPDF2.PdfWriter()

        # Extract the first pages of each PDF
        base_page = base_reader.pages[0]
        overlay_page = overlay_reader.pages[0]

        # Merge overlay onto the base page
        base_page.merge_page(overlay_page)

        # Add the merged page to the writer object
        writer.add_page(base_page)

        # Save the final merged PDF
        with open(output_pdf, "wb") as out_file:
            writer.write(out_file)


# Main function managing the pipeline execution.
# Parses command-line arguments, processes BAM and event files, generates Sashimi plots, 
# and integrates annotations for visualization.
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-directory', dest='directory', help='Path to output directory')
    parser.add_argument('-ggsashimi', dest='ggsashimi', help='Full path to executable ggsashimi')
    parser.add_argument('-bam', dest='bam', help='Full path to BAM file of interest')
    parser.add_argument('-event_file', dest='event_file', help='Full path to event file of interest')
    parser.add_argument('-extend_bp', dest='extend_bp', help='Number of bp to extend at extremities. Default = 50', type=int, default=50)
    parser.add_argument('-MinThresholdNbReads', dest='MinThresholdNbReads', help='Minimum reads in junctions for plot. Default: Half of the reads in the sample of interest junction', type=int, default=-1)
    parser.add_argument('-list_bam', dest='list_bam', help='List of BAM files')
    parser.add_argument('-nb_samples', dest='nb_samples', help='Number of random samples in plot. Default = 4', type=int, default=4)
    parser.add_argument('-gtf', dest='gtf', help='Full path to reference GTF file', default=-1)
    parser.add_argument('-color', dest='color', help='Palette of colors for plot', default=-1)
    args = parser.parse_args()

    # Extract input parameters
    directory = args.directory
    ggsashimi = args.ggsashimi
    bam = args.bam
    name_bam = bam.rsplit("/",1)[1].rsplit(".",2)[0]
    event_file = args.event_file
    extend_bp = args.extend_bp
    MinThresholdNbReads = args.MinThresholdNbReads
    list_bam = args.list_bam.split() if args.list_bam else []
    nb_samples = args.nb_samples
    gtf = args.gtf
    color = args.color
    python = sys.executable

    # Ensure the output directory exists
    create_directory_if_not_exists(directory)

    # Read event file (CSV or Excel format)
    extension = event_file.rsplit(".", 1)[1]
    df = pd.read_excel(event_file, engine='openpyxl') if extension == "xlsx" else pd.read_csv(event_file, sep="\t")
    if df.empty:
        print(f"No events found. The file {event_file} is empty.")
        return

    # Display execution summary
    print("######################")
    print(f'# Generating Sashimi plot for {name_bam}...')
    print("######################")
    
    # Iterate through event entries
    for index, row in df.iterrows():
        directory = args.directory
        chr = row["chr"]
        start = int(row["start"])
        end = int(row["end"])
        gene = str(row["Gene"])
        event_type = str(row["event_type"])
        nm = str(row["NM"])
        filterInterpretation = str(row["filterInterpretation"])

        if filterInterpretation == "No model" :
            directory = directory + "/no_model"
            create_directory_if_not_exists(directory)
        
        if filterInterpretation == "Percentage threshold execeeded" :
            directory = directory + "/threshold_exceeded"
            create_directory_if_not_exists(directory)

        if filterInterpretation == "Unique junction" :
            directory = directory + "/unique"
            create_directory_if_not_exists(directory)

        if filterInterpretation == "Event too complex" :
            directory = directory + "/event_too_complex"
            create_directory_if_not_exists(directory)

        if gtf != -1:
            pdf_temp = directory + "/" + row["Conca"] + "_temp.pdf"
        else:
            pdf_temp = directory + "/" + row["Conca"] + ".pdf"

        # Determine threshold for reads
        if MinThresholdNbReads == -1:
            MinThresholdNbReads = int(row[[name_bam]].iloc[0] / 2)

        # # Select random BAM files for visualization
        bams = ""

        # Ajout du BAM principal
        bam = os.path.abspath(bam)
        bams += f"{name_bam}\t{bam}\tInterest\n"

        # Préparation de la liste sans le BAM principal
        remaining_bams = []
        for path in list_bam:
            name = path.rsplit("/", 1)[1].rsplit(".", 2)[0]
            if name != name_bam:
                remaining_bams.append(path)

        # Tirage aléatoire explicite
        for i in range(nb_samples):
            if not remaining_bams:
                break  # sécurité si la liste est trop courte
            index = random.randint(0, len(remaining_bams) - 1)
            selected = os.path.abspath(remaining_bams.pop(index))
            name = selected.rsplit("/", 1)[1].rsplit(".", 2)[0]
            bams += f"{name}\t{selected}\tRandom\n"
        if bams.endswith("\n"):
            bams = bams[:-1]

        success = False
        # Tentatives avec élargissement progressif : 1x, 2x, 3x, 4x, 5x, 6x
        for attempt, factor in enumerate([1, 2, 3, 4, 5, 6], start=1):
            start_extend = start - extend_bp * factor
            end_extend   = end + extend_bp * factor

            cmd = f""" bash -c '{python} {ggsashimi} -b <(echo -e "{bams}") \
                -c {chr}:{start_extend}-{end_extend} \
                -M {MinThresholdNbReads} --fix-y-scale --height 2.75 --ann-height 2.75 \
                --alpha 0.5 -o {pdf_temp} """

            if str(gtf) != "-1":
                cmd += f" -g {gtf}"
            if str(color) != "-1":
                cmd += f" -P {color} -C 3"
            cmd += "'"
            MinThresholdNbReads = args.MinThresholdNbReads
            print(f"[INFO] Tentative {attempt} avec fenêtre ±{extend_bp*factor} bp")
            ret = os.system(cmd)

            if ret == 0 and os.path.exists(pdf_temp):
                print(f"[SUCCESS] Sashimi plot généré : {pdf_temp}")
                success = True
                break
            else:
                print(f"[WARN] Échec tentative {attempt} pour {row['Conca']}")

        if not success:
            print(f"[ERROR] Impossible de générer le sashimi plot pour {row['Conca']} après 6 tentatives")
            continue



        # Annotate Sashimi plots with transcript labels if GTF is provided
        if gtf != -1:
            pdf_final = directory + "/" + row["Conca"] + ".pdf"
            overlay_pdf = directory + "/" + row["Conca"] + ".overlay.pdf"
            region = f"{chr}:{start_extend}-{end_extend}"
            
            expanded_pdf = pdf_temp.replace(".pdf", "_expanded.pdf")
            N = 3 #Extend


            with open(pdf_temp, "rb") as original_file, open(expanded_pdf, "wb") as new_file:
                reader = PyPDF2.PdfReader(original_file)
                writer = PyPDF2.PdfWriter()

                for page in reader.pages:
                    # Extend PDF to the left
                    page.mediabox.lower_left = (page.mediabox.lower_left[0] - N, page.mediabox.lower_left[1])
                    page.mediabox.upper_right = (page.mediabox.upper_right[0], page.mediabox.upper_right[1]) 

                    # Adjust CropBox to display space correctly
                    page.cropbox.lower_left = (page.cropbox.lower_left[0] - N, page.cropbox.lower_left[1])
                    page.cropbox.upper_right = (page.cropbox.upper_right[0], page.cropbox.upper_right[1])

                    writer.add_page(page)
                
                writer.write(new_file)

            os.remove(pdf_temp)
            pdf_temp = expanded_pdf

            exon_annotations = read_gtf_and_number_exons(gtf, region)
            pdf_width, pdf_height = get_pdf_dimensions(pdf_temp)
            positions, max_y = calculate_transcript_labels(pdf_temp, exon_annotations, pdf_width, pdf_height, nm, (1 -N), gtf)
            
            create_labels_overlay(overlay_pdf, positions, pdf_width, pdf_height, chr, start, end, gene, event_type, max_y, (1 -N))
            reader_overlay = PyPDF2.PdfReader(overlay_pdf)
            writer = PyPDF2.PdfWriter()

            for page in reader_overlay.pages:
                # Move overlay repository to -N
                page.mediabox.lower_left = (page.mediabox.lower_left[0] - N, page.mediabox.lower_left[1])
                page.mediabox.upper_right = (page.mediabox.upper_right[0] - N, page.mediabox.upper_right[1])

                writer.add_page(page)

            overlay_adjusted = overlay_pdf.replace(".pdf", "_adjusted.pdf")
            writer.write(overlay_adjusted)
            merge_pdfs(pdf_temp, overlay_adjusted, pdf_final)
            
            print(f"Annotated PDF created: {pdf_final}")
            
            os.remove(pdf_temp)
            os.remove(overlay_adjusted)
            os.remove(overlay_pdf)

# Entry point for script execution
if __name__ == "__main__":
    main()