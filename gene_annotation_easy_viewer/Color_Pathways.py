import pickle
import os
import sys
import colorsys
import math
import copy
import tkinter as tk
from collections import namedtuple
from tkinter import filedialog, messagebox, simpledialog
from tkinter import *
from GAEV import *

cur_dir = os.getcwd()

sets_to_highlight = {}  # {"blue": [gene1, gene,2 gene5], "green": [gene4, gene7, gene9], ...}
pathway_info_full = {}
pathway_info_query = {}
gene_info_full = {}
gene_info_query = {}
color_hex_dict = {}  # will store the color hex of each k code in query {"K14515": }


def load_full_annotation(data_file_path_full):
    # imports data from files and stores it into appropriate lists and dict for later access
    with open(data_file_path_full, "rb") as f:  # open data file that was generated previously
        if not pickle.load(f):  # reads the first line of saved data which tells whether it is complete or not
            messagebox.showerror("Data is not complete")  # exits program if data is not complete
            return

        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in genes
            gene = pickle.load(f)  # stores the Gene object as gene
            if gene.k_code in gene_info_full:  # dictionary with the k code as the key and gene object stored in list
                gene_info_full[gene.k_code].append(gene)  # adds gene to list of other genes with shared k codes in dict
            else:
                gene_info_full[gene.k_code] = [gene]  # adds new k code entry in dict as key with corresponding gene in list
        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in pathways
            pathway = pickle.load(f)
            pathway_info_full[pathway.map_code] = pathway  # loads data into pathway_info_1


# will simply read file and return list with each striped line on an element
def load_set(set_input_path):
    output_list = []
    with open(set_input_path, 'r') as f:
        for line in f:
            output_list.append(line.strip())  # strip to remove unexpected blank spaces / new lines

    return output_list


def generate_gene_info_query():  # generates a deep copy of gene_info_full and only keeps genes that are in sets_to_highlight
    global gene_info_query

    # makes a new copy of gene_info_full that can be changed independently of the original
    gene_info_query = copy.deepcopy(gene_info_full)
    # removes all genes in gene_info_query that are not in the set to be highlighted, removed genes will appear as
    # grey on the pathway
    for k_code in gene_info_full:
        gene_info_query[k_code] = [gene for gene in gene_info_full[k_code]
                                   if any([gene.gene_num in gene_set for gene_set in sets_to_highlight.values()])]

    # removes all k_codes form gene_info_query dict that has no genes in the set to highlight
    gene_info_query = dict([(k_code, gene_list) for k_code, gene_list in gene_info_query.items() if gene_list])


def generate_pathway_info_query():
    global pathway_info_query

    # makes a new copy of pathway_info_full that can be changed independently of the original
    pathway_info_query = copy.deepcopy(pathway_info_full)

    # empties list of genes involved in pathway from every pathway in query list, so that it may be repopulated with
    # only genes from gene_info_query
    for map_code in pathway_info_full:
        pathway_info_query[map_code].genes_invol = []

    for k_code in gene_info_query:
        # only have to iterate through one gene object for each k_code since genes with the same k_code will have
        # identical  link_path
        for map_info in gene_info_query[k_code][0].link_path:  # map info is [m-code\description, ...]
            m_code = map_info[:8]    # isolates the map code which is always 8 characters map#####
            # calls method to add k_code to genes_invol while keeping the k_code ordered and preventing duplicates
            pathway_info_query[m_code].add_gene(k_code)

    pathway_info_query = dict([(m_code, pathway) for m_code, pathway in pathway_info_query.items()
                               if pathway.genes_invol and len(pathway_info_full[m_code].genes_invol) < 104])


# will specific the color of each gene by K code, if more than one color is used then
def generate_color_hex_dict():
    for k_code in gene_info_query:
        color_list = []
        # finds highlight color associated with the gene and adds it to the color list
        for gene in gene_info_query[k_code]:
            for color in sets_to_highlight.keys():
                if gene.gene_num in sets_to_highlight[color]:
                    color_list.append(color)

        # will blend all colors from all genes associated with the k_code
        # ex. blue + red = purple; it is weighted by occurrence of color blue + 3 red = dark pink
        blend_dict = {}  # {color_hex: weight, color_hex: weight, ...} where weight is number of occurrence
        for color in set(color_list):
            blend_dict[color] = color_list.count(color)
        color_hex_dict[k_code] = combine_hex_values(blend_dict)

    # for every k_code that is in the full annotation, but not included in any sets to highlight
    for k_code in gene_info_full:
        if k_code not in color_hex_dict.keys():
            color_hex_dict[k_code] = "d9d9d9"  # stores a hex for a light grey color


# accepts a dict of {color_hex: weight, color_hex: weight, ...} and returns a color hex that blends all the colors
def combine_hex_values(d):
    d_items = sorted(d.items())
    tot_weight = sum(d.values())
    red = int(sum([int(k[:2], 16)*v for k, v in d_items])/tot_weight)
    green = int(sum([int(k[2:4], 16)*v for k, v in d_items])/tot_weight)
    blue = int(sum([int(k[4:6], 16)*v for k, v in d_items])/tot_weight)
    zpad = lambda x: x if len(x)==2 else '0' + x
    return zpad(hex(red)[2:]) + zpad(hex(green)[2:]) + zpad(hex(blue)[2:])


def _new_generate_url(self):  #
    self.url = "http://www.kegg.jp/kegg-bin/show_pathway?map=" + self.map_code + "&multi_query="  # sets the which pathway to use

    for k_code in self.genes_invol:  # goes through each unique k_code in pathway
        self.url = self.url + k_code + "+%23" + color_hex_dict[k_code] + "%0a"  # adds gene's k_code to the end of the url and specifies color

    return self.url
Pathway_MAP.generate_url = _new_generate_url  # overrides the generate_url method to include unique color value


def make_html_table_rows(sorted_m_codes):
    html_rows = ""
    for m_code in sorted_m_codes:
        # generates the text for the hyperlink (name + num of associated genes + percentage of genes out of total)
        query_pathway = pathway_info_query[m_code]
        full_pathway = pathway_info_full[m_code]
        hyper_text = query_pathway.name + "(" + str(len(query_pathway.genes_invol)) + ", " + \
                     "{:.2f}%".format(len(query_pathway.genes_invol) / len(full_pathway.genes_invol) * 100) + ")"
        url = full_pathway.generate_url()  # generates the url
        hyperlink = "<a href=\"" + url + "\">" + hyper_text + "</a>"  # embeds hyperlink to text
        row = "<tr><td>" + hyperlink + "</td></tr>"
        html_rows = html_rows + row
    return html_rows


class PathwayInfo:
    query_pathway: Pathway_MAP
    full_pathway: Pathway_MAP
    percentage: float

    def __init__(self, query_pathway, full_pathway, color_hex, percentage=0):
        self.query_pathway = query_pathway
        self.full_pathway = full_pathway
        self.percentage = percentage


def generate_html(output_path):
    with open(output_path, 'w+') as f:
        html_code = """<html>
                <head>
                <style>
                    html {
                        /*width: 100%;*/ /*required if using % width*/
                        /*height: 100%;*/ /*required if using % height*/
                    }
                    body {
                        /*width: 100%;*/ /*required if using % width*/
                        /*height: 100%;*/ /*required if using % height*/
                        /*margin: 0;*/ /*required if using % width or height*/
                        /*padding: 0 20px 0 20px;*/ /*purely aesthetic, not required*/
                      /*box-sizing: border-box;*/ /*required if using above declaration*/
                      background: white;
                        text-align: center; /*delete if using % width, or you don't want table centered*/
                    }
                    .scrollingtable {
                        box-sizing: border-box;
                        display: inline-block;
                        vertical-align: middle;
                        overflow: hidden;
                        width: auto; /*set table width here if using fixed value*/
                        /*min-width: 100%;*/ /*set table width here if using %*/
                        height: 600px; /*set table height here; can be fixed value or %*/
                        /*min-height: 104px;*/ /*if using % height, make this at least large enough to fit scrollbar arrows + captions + thead*/
                        font-family: Verdana, Tahoma, sans-serif;
                        font-size: 15px;
                        line-height: 20px;
                        padding-top: 20px; /*this determines top caption height*/
                        padding-bottom: 20px; /*this determines bottom caption height*/
                        text-align: left;
                    }
                    .scrollingtable * {box-sizing: border-box;}
                    .scrollingtable > div {
                        position: relative;
                        border-top: 1px solid black; /*top table border*/
                        height: 100%;
                        padding-top: 20px; /*this determines column header height*/
                    }
                    .scrollingtable > div:before {
                        top: 0;
                        background: cornflowerblue; /*column header background color*/
                    }
                    .scrollingtable > div:before,
                    .scrollingtable > div > div:after {
                        content: "";
                        position: absolute;
                        z-index: -1;
                        width: 100%;
                        height: 50%;
                        left: 0;
                    }
                    .scrollingtable > div > div {
                        /*min-height: 43px;*/ /*if using % height, make this at least large enough to fit scrollbar arrows*/
                        max-height: 100%;
                        overflow: scroll; /*set to auto if using fixed or % width; else scroll*/
                        overflow-x: hidden;
                        border: 1px solid black; /*border around table body*/
                    }
                    .scrollingtable > div > div:after {background: white;} /*match page background color*/
                    .scrollingtable > div > div > table {
                        width: 100%;
                        border-spacing: 0;
                        margin-top: -20px; /*inverse of column header height*/
                        /*margin-right: 17px;*/ /*uncomment if using % width*/
                    }
                    .scrollingtable > div > div > table > caption {
                        position: absolute;
                        top: -20px; /*inverse of caption height*/
                        margin-top: -1px; /*inverse of border-width*/
                        width: 100%;
                        font-weight: bold;
                        text-align: center;
                    }
                    .scrollingtable > div > div > table > * > tr > * {padding: 0;}
                    .scrollingtable > div > div > table > thead {
                        vertical-align: bottom;
                        white-space: nowrap;
                        text-align: center;
                    }
                    .scrollingtable > div > div > table > thead > tr > * > div {
                        display: inline-block;
                        padding: 0 6px 0 6px; /*header cell padding*/
                    }
                    .scrollingtable > div > div > table > thead > tr > :first-child:before {
                        content: "";
                        position: absolute;
                        top: 0;
                        left: 0;
                        height: 20px; /*match column header height*/
                        border-left: 1px solid black; /*leftmost header border*/
                    }
                    .scrollingtable > div > div > table > thead > tr > * > div[label]:before,
                    .scrollingtable > div > div > table > thead > tr > * > div > div:first-child,
                    .scrollingtable > div > div > table > thead > tr > * + :before {
                        position: absolute;
                        top: 0;
                        white-space: pre-wrap;
                        color: white; /*header row font color*/
                    }
                    .scrollingtable > div > div > table > thead > tr > * > div[label]:before,
                    .scrollingtable > div > div > table > thead > tr > * > div[label]:after {content: attr(label);}
                    .scrollingtable > div > div > table > thead > tr > * + :before {
                        content: "";
                        display: block;
                        min-height: 20px; /*match column header height*/
                        padding-top: 1px;
                        border-left: 1px solid black; /*borders between header cells*/
                    }
                    .scrollingtable .scrollbarhead {float: right;}
                    .scrollingtable .scrollbarhead:before {
                        position: absolute;
                        width: 100px;
                        top: -1px; /*inverse border-width*/
                        background: white; /*match page background color*/
                    }
                    .scrollingtable > div > div > table > tbody > tr:after {
                        content: "";
                        display: table-cell;
                        position: relative;
                        padding: 0;
                        border-top: 1px solid black;
                        top: -1px; /*inverse of border width*/
                    }
                    .scrollingtable > div > div > table > tbody {vertical-align: top;}
                    .scrollingtable > div > div > table > tbody > tr {background: white;}
                    .scrollingtable > div > div > table > tbody > tr > * {
                        border-bottom: 1px solid black;
                        padding: 0 6px 0 6px;
                        height: 20px; /*match column header height*/
                    }
                    .scrollingtable > div > div > table > tbody:last-of-type > tr:last-child > * {border-bottom: none;}
                    .scrollingtable > div > div > table > tbody > tr:nth-child(even) {background: gainsboro;} /*alternate row color*/
                    .scrollingtable > div > div > table > tbody > tr > * + * {border-left: 1px solid black;} /*borders between body cells*/
                </style>
                </head>
                <body>
                <body link="#003366">"""

        html_code = html_code + """
                <br />
                <div class="scrollingtable">
                    <div>
                        <div>
                            <table>
                                <caption>Pathway Sorted by Number of Associated Genes</caption>
                                <thead>
                                    <tr>
                                        <th>
                                            <!--more versatile way of doing column label; requires 2 identical copies of label-->
                                            <div><div>Pathway</div><div>Pathway</div></div>
                                        </th>
                          <th class="scrollbarhead"></th> <!--ALWAYS ADD THIS EXTRA CELL AT END OF HEADER ROW-->
                                    </tr>
                                </thead>
                                <tbody>
                                    """
        # sorts m_code by the number of highlighted genes in the table to iterate over dictionary in desired order
        # determines the order of the pathways in the html table
        sorted_m_code = sorted(pathway_info_query.keys(),
                               key=lambda m: len(pathway_info_query[m].genes_invol), reverse=True)
        html_code = html_code + make_html_table_rows(sorted_m_code)
        html_code = html_code + """</tbody>
                            </table>
                        </div>
                    </div>
                </div>
                <span style="display:inline-block; width: 50;"></span>
                <div class="scrollingtable">
                    <div>
                        <div>
                            <table>
                                <caption>Pathway Sorted by Name</caption>
                                <thead>
                                    <tr>
                                        <th>
                                            <!--more versatile way of doing column label; requires 2 identical copies of label-->
                                            <div><div>Pathway</div><div>Pathway</div></div>
                                        </th>
                          <th class="scrollbarhead"></th> <!--ALWAYS ADD THIS EXTRA CELL AT END OF HEADER ROW-->
                                    </tr>
                                </thead>
                                <tbody>
                                    """
        # sorts m_code alphabetically by pathway name to iterate over dictionary in desired order; determines the order
        # of the pathways in the html table
        sorted_m_code = sorted(pathway_info_query.keys(),
                               key=lambda m: pathway_info_query[m].name)
        html_code = html_code + make_html_table_rows(sorted_m_code)
        html_code = html_code + """</tbody>
                            </table>
                        </div>
                    </div>
                </div>
                </body>
                </html>"""
        f.write(html_code)
    f.close()


class UI:
    color_dict: dict
    fields_dict: dict

    def __init__(self, root):
        # sets the color palette by creating dict with name as key and hex value as value
        self.color_dict = {"blue": "7DF9FF", "purple": "DB70FF", "red": "FF261F", "orange": "FFAE42", "yellow": "FFF600",
                           "green": "39FF14", "pink": "FF66CC", "brown": "AF6E4D"}

        self.fields_dict = {}
        self.root = root
        self.root.title("Highlight Pathways")
        self.root.geometry('800x300')
        self.instantiate_window()
        self.root.mainloop()

    def text_entry(self, field):
        Label(self.root,
              text="Enter path to the GAEV annotation reference file. This is the dat file produced by GAEV from the full KAAS annotation list."
              ).pack(side=TOP, anchor=W, padx=5, pady=5)

        row = Frame(self.root)

        self.fields_dict[field] = StringVar(self.root)
        button = Button(row, text="Browse...", command=lambda: self.set_reference_annotation(field))

        label = Label(row, width=20, text=field, anchor='w')
        entry = Entry(row, textvariable=self.fields_dict[field])
        row.pack(side=TOP, fill=X, padx=5, pady=5)
        label.pack(side=LEFT)
        entry.pack(side=LEFT, expand=YES, fill=X)
        button.pack(side=LEFT, padx=2)

    def instantiate_window(self):
        self.text_entry("Full GEAV Annotation File")

        row = Frame(self.root)
        add_set_frame = Frame(row)
        Button(add_set_frame, text="Add a new gene set to highlight",
               command=lambda: self.input_new_set()).pack(side=RIGHT, padx=2)
        add_set_frame.pack(side=LEFT, expand=True, fill=X)
        export_frame = Frame(row)
        Button(export_frame, text="Export pathways",
               command=lambda: self.output_menu()).pack(side=LEFT, padx=2)
        export_frame.pack(side=LEFT, expand=True, fill=X)
        row.pack(side=TOP, fill=X, padx=5, pady=(20,5))

        curr_set_table = Frame(self.root)
        self.fields_dict["sets"] = StringVar(self.root)
        set_frame = Frame(curr_set_table)
        set_frame.pack(side=LEFT, expand=True, padx=15, fill=BOTH)
        Label(set_frame, text="Current sets:").pack(side=TOP, anchor=W, padx=5)
        Label(set_frame, textvariable=self.fields_dict["sets"], bg="white", anchor=N).pack(side=LEFT, expand=True,
                                                                                                padx=5, pady=5, fill=BOTH)
        color_frame = Frame(curr_set_table)
        color_frame.pack(side=LEFT, expand=True, padx=15, fill=BOTH)
        self.fields_dict["color"] = StringVar(self.root)
        Label(color_frame, text="Color:").pack(side=TOP, anchor=W, padx=5)
        Label(color_frame, textvariable=self.fields_dict["color"], bg="white", anchor=N).pack(side=LEFT, expand=True,
                                                                                                 padx=5, pady=5, fill=BOTH)
        curr_set_table.pack(side=TOP, expand=True, fill=BOTH, pady=(5,15))


    def set_reference_annotation(self, field):
        data_file_path_full = filedialog.askopenfilename(initialdir=os.getcwd(),
                                                title="Enter path to annotation reference file. This is the dat file produced by GAEV on the full annotation list.",
                                                filetypes=(("GAEV output files", "*.dat*"), ("all files", "*.*")))
        os.chdir(cur_dir)

        try:
            load_full_annotation(data_file_path_full)
            self.fields_dict[field].set(data_file_path_full)
        except:
            messagebox.showerror(title="Error", message="Could not process GAEV dat file.")

    def input_new_set(self):
        if not gene_info_full:
            messagebox.showwarning(title="Warning", message="Please load a GAEV annotation reference file before proceeding.")
            return

        input_path = filedialog.askopenfilename(initialdir=os.getcwd(),
                                                title="Please specify path to gene set. Input file should contain one gene ID per line.")
        if not input_path:
            return

        color_hex = self.input_color()
        if not color_hex:
            return

        self.fields_dict["sets"].set(self.fields_dict["sets"].get() + input_path + "\n")

        # adds list of gene ids and the color of those genes to the global list
        sets_to_highlight[color_hex] = load_set(input_path)

    def input_color(self):  # prompts for color to set the fill color of the genes in the set on the pathway maps
        color_text = simpledialog.askstring(title="Choose color",
                                          prompt="Enter the highlight color for this set. Use the name of the color or input the color hex")

        if color_text.lower() in self.color_dict:  # will accept color names (case-insensitive) set in color_dict
            self.fields_dict["color"].set(self.fields_dict["color"].get() + color_text + "\n")
            return self.color_dict[color_text.lower()]
        elif len(color_text) == 6:  # else it will check if the input was a color hex
            try:
                int(color_text.upper(), 16)
                self.fields_dict["color"].set(self.fields_dict["color"].get() + color_text + "\n")
                return color_text.upper()
            except ValueError:
                pass

        print(textwrap.dedent("""
                                 Color entered is not valid!"""))  # displays options for
        return self.input_color()

    def output_menu(self):

        if not gene_info_full:
            messagebox.showwarning(title="Warning", message="Please load a GAEV annotation reference file before proceeding.")
            return

        output_path = filedialog.asksaveasfilename(initialdir=cur_dir,
                                                   filetypes=[('html', '*.html'), ('All Files', '*.*')],
                                                   defaultextension='.html')
        print(output_path + '\n')

        generate_gene_info_query()
        generate_pathway_info_query()
        generate_color_hex_dict()
        generate_html(output_path)

        self.root.destroy()


if __name__ == "__main__":
    root = tk.Tk()
    ui = UI(root)
