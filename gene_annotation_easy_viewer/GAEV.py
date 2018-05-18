############ Imports ############
import contextlib  # the closing context manager will ensure connection to url is closed after with statement
import textwrap  # used when printing UI in order to remove common leading indents
import os  # used to remove temporary files created and path manipulation
import urllib  # to handle HTTPError exceptions raised from the urllib module
from urllib.request import urlopen  # import the tools we need to open url
import bisect  # module needed to insert an item into a sorted list
import pickle  # needed to save and load data
import sys  # needed to exit out of program when error is encountered
import re  # required to split files without removing delimiter

############ Variables ############
_input_file = "No_File_Specified"  # will store path of the input file
_trimmed_file = "No_File_Specified"  # will store path of the trimmed input file with genes not associated with KO# removed
_data_file = "No_File_Specified" # will store path to the data file code will generate
_html_file = "No_File_Specified" # will store path to the html file that will display the data
_pathway_list = []  # will store pathway list loaded from data file
_gene_list = []  # will store gene list loaded from data file
_total_genes = 0  # will store the total number of genes when unfiltered

############ Classes and Functions ############
def decode_url(urlLink):  # converts HTML response into String (allows program to read webpages)
    while True:  # will loop forever until connection with url is made
        with contextlib.closing(urlopen(urlLink)) as response:  # opening the url gives us HTML response variables instead of String
            html_response = response.read()
            encoding = response.headers.get_content_charset('utf-8')  # handles the encoding from Content-Type
            decoded_html = html_response.decode(encoding)
        response.close()
        return decoded_html  # returns the decoded html as String

def save(geneList, pathwayList, completed):
    global _data_file  # sets local _data_file to global _data_file
    with open(_data_file, "wb") as f:  # create a file to save/pickle data
        pickle.dump(completed, f)
        pickle.dump(len(geneList), f)  # store the amount of data entries for genes
        for gene in geneList:  # loops through and store all all data in geneList
            pickle.dump(gene, f)
        pickle.dump(len(pathwayList), f)  # store the amount of data entries for pathways
        for pathway in pathwayList:  # loops through and store all data in pathwayList
            pickle.dump(pathway, f)
    f.close()

def separate_file(path):  # accepts path as string and separates the file name ['C:/ExampleFolder/', 'FileName']
    path_list = re.split("(/)", path)  # separates the string at every "/" while keeping the delimiter
    return [path_list[:-1],path_list[-1]]  # returns list with all the path without file name and file name

def get_file_name(path):  # returns the file new from input path (should work on all operating system)
    head, tail = os.path.split(path)  # splits path into directory (head) and file name (tail)
    file_name = tail or os.path.basename(head)  # if path ended in slash (a\b\c\) will use 'basename()' to get file name
    file_no_ext, ext = os.path.splitext(file_name)  # will separate ext from file name
    return file_name, file_no_ext  # will return a tuple of file name with ext followed by file name without ext

def set_input(path):  # used to set the paths of input file, trimmed file, and data file based on input file path
    global _input_file  # sets local _input_file to global _input_file
    global _trimmed_file  # sets local _trimmed_file to global _trimmed_file
    global _data_file  # sets local _data_file to global _data_file
    _input_file = os.path.abspath(path)  # handles distinction between relative and abs path by converting all to abs path
    file_dir = os.path.dirname(_input_file)
    file_name, file_no_ext = get_file_name(_input_file)  # gets file name (w/ and w/o ext) without the directory  path
    _trimmed_file = os.path.join(file_dir, "trimmed_" + file_name)  # create path for trimmed file with dir and file name
    _data_file = os.path.join(file_dir, file_no_ext + ".dat")  # create path for data file with ".dat" extension
    set_data_file()  # with no argument passed, functions will just store a path to html using the stored _data_file


def set_data_file(data_file = None):  # will specify path to data file to be used and path of output html file
    global _data_file  # sets local _data_file to global _data_file
    global _html_file  # sets local _html_file to global _html_file

    if data_file == None:  # if a path to the data_file was not specify,
        data_file = _data_file  # use the path that was most recently stored
    else:  # if a path was specified,
        data_file = os.path.abspath(data_file)  # obtains the absolute path to the specified data file
        _data_file = data_file  # then store it as the global _data_file path
        global _trimmed_file  # sets local _trimmed_file to global _trimmed_file
        global _input_file  # sets local _input_file to global _input_file
        _input_file = "No_File_Specified"  # (1/2)  delete the global _input_file and the _trimmed_file path to avoid
        _trimmed_file = "No_File_Specified"  # (2/2) overriding the data file with the wrong dataset

    name_no_ext, ext = os.path.splitext(data_file)  # stores file name w/o ext in name_no_ext and the extension in ext
    _html_file = os.path.join(os.path.dirname(data_file),  name_no_ext + ".html")  # creates path to html file

def load_data(data_file = None):  # will load data from data file into _gene_list and _pathway_list
    global _data_file  # sets local _data_file to global _data_file
    global _pathway_list  # sets local _pathway_list to global _pathway_list
    global _gene_list  # sets local _gene_list to global _gene_list
    global _total_genes  # sets local _total_genes to global _gene_list

    del _pathway_list[:]  # clears _pathway_list to prevent appending one data set onto another
    del _gene_list[:]  # clears _gene_list to prevent appending one data set onto another

    if data_file == None:  # if a path to the data_file was not specify,
        data_file = _data_file # use the path that was most recently stored
    with open(data_file, "rb") as f:  # open data file that was generated previously
        if not pickle.load(f):  # reads the first line of saved data which tells whether it is complete or not
            sys.exit("Data is not complete")  # exits program if data is not complete

        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in genes
            _gene_list.append(pickle.load(f))  # loads data into _gene_list
        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in pathways
            _pathway_list.append(pickle.load(f))  # loads data into _pathway_list

    _total_genes = len(_gene_list)

def retrieve_batch_files(data):
    file_list = []  # will be used to store list of input/data files
    for file_path in os.listdir('.'):  # cycle through every file in the current directory, but...
        if data:  # if user specified they wanted to use data files
            if file_path.endswith(".dat"):  # only use files with the .dat extension
                file_list.append(file_path)
        else:  # if user specified they wanted to use input files
            # only use txt files, but do not use any of the output txt files as a precaution
            if file_path.endswith(".txt") and not file_path.endswith("-gene_table.txt"):
                file_list.append(file_path)  # add the input file path to the input_file_list
    return file_list

class Pathway_MAP:  # class for pathway map objects
    def __init__(self, ipathway_info):  # accepts the map code of the pathway
        pathway_info = ipathway_info.split(None, 1)  # separtes the map code and the map name by spiting at first space
        self.map_code = pathway_info[0]  # stores map code in map_code
        self.name = pathway_info[1]  # stores name into the 'name' variable
        self.genes_invol = []  # creates an empty list that will store all the genes involved in the pathway
        self.url = ""  # base url or pathway map with genes highlighted

    def add_gene(self, ik_code):
        bisect.insort(self.genes_invol, ik_code)  # add a gene to pathway and keeps gene ordered by k_code

    def generate_url(self, ik_code):  # generates base url from ko code of pathway and k code of genes involved and selected gene
        self.url = "http://www.kegg.jp/kegg-bin/show_pathway?map=" + self.map_code + "&multi_query="  # sets the which pathway to use
        for k_code in self.genes_invol:  # goes through each gene in genes_invol
            self.url = self.url + k_code + "+%23bfffbf%0d%0a"  # adds gene's k_code to the end of the url and specifies color
        self.url = self.url + ik_code + "+%238B0000,%23F0F8FF"  # adds the selected gene and marks it with a different color
        return self.url

    def check_genes(self, gene_list):  # checks if any associated genes is in the current filtered gene list
        k_code = [gene.k_code for gene in gene_list]  # extracts k-codes from list of gene objects into their own list
        is_there = bool(set(self.genes_invol) & set(k_code))  # tests for common items between genes_invol and k_code
        return is_there


class Gene:  # class for gene objects
    def __init__(self,ig_num, ik_code):  # accepts the gene number and the KEGG's k code for the gene
        self.gene_num = ig_num  # sets gene number
        self.k_code = ik_code  # sets k number
        infoList = self.get_info()  # runs get_info() to retrieve information on gene from KEGG
        self.name = infoList[0]  # sets gene name
        self.definition = infoList[1]  # sets gene definition
        self.link_path = infoList[2]  # sets pathways that gene is involved in [map_code, map_name, map_code, map_name ...]

    def format_pathway_info(self, ipathway_info):  # accepts line from gene's "get" site from KEGG
        ipathway_info = "map" + ipathway_info[2:]  # replaces the 'ko' prefix with 'map' prefix
        return ipathway_info

    def get_info(self):  # accepts a k code and returns a list with name first then description then associated map codes
        infoList = []  # empty list to store gene name and description and pathways the gene is involved in
        # ex list: [NAME, DEFINITION, [map_code map_name, map_code map_name, map_code map_name, ...]]
        linkPathList = []  # list to store pathways the gene is involved in [map_code map_name, ...]
        link = 'http://rest.kegg.jp/get/' + self.k_code  # creates link to get info on gene
        try:
            info_url = decode_url(link)  # access KEGG for information
        except urllib.error.HTTPError as error:
            return ['', '', []]
        for line in info_url.splitlines():  # splits string by line
            if not linkPathList:  # pythonic way to check if code reached pathway yet, if it has then go to else
                column = line.split(None, 1)  # split line only once at first space to create two columns
                if column[0] == 'NAME':  # checks first column to see if it is a name by checking first column
                    # ex line: NAME ITPR1
                    infoList.append(column[1])  # adds second column (actual name) to infoList; always first
                if column[0] == 'DEFINITION':  # checks first column to see if it is a definition by checking first column
                    # ex line: DEFINITION inositol 1,4,5-tripohsphate receptor type 1
                    infoList.append(column[1])  # adds second column (definition) to infoList; always second
                if column[0] == 'PATHWAY':  # checks first column to see if list of pathways is starting by checking first column
                    # ex line: Pathway ko04020 Calcium signaling pathway
                    linkPathList.append(self.format_pathway_info(column[1]))  # adds formatted pathway info to linkPathList
            else:  # if started adding linked pathways, just keep with the chug and plug
                line = line.lstrip()  # strips all the blank spaces from the beginning of the line
                if line[:2] == "ko":  # checks if first two characters is ko to comfirm still on pathway
                    linkPathList.append(self.format_pathway_info(line))  # adds formatted pathway info in linkPathList
                    # ex line: ko00010  Glycolysis / Gluconeogenesis
                else:  # after all pathways are added
                    break  # break out of for loop because there is no reason to keep checking file
        infoList.append(linkPathList)  # add linkPathList to infoList
        return infoList  # finally return infoList [NAME, DEFINITION, [map_code map_name, map_code map_name, ...]]

    # searches name to check if target term is in the name
    def check_name(self, target):
        is_there = False  # boolean to check if target term is in name
        name = self.name
        if target.lower() in name.lower():  # checks if name contains target term
            is_there = True
        return is_there

    # searches definition to check if target term is in the definition
    def check_definition(self, target):
        is_there = False  # boolean to check if target term is in definition
        definition = self.definition
        if target.lower() in definition.lower():  # checks if definition contains target term
            is_there = True
        return is_there

    #  searches linked pathways to see if it includes the target term
    def check_pathway(self, target, pathway_list):
        is_there = False  # boolean to check if target term is in linked pathways
        for m_code in self.link_path:  # cycles through every linked pathway in gene
            # each element in link_path is ['map_code map_name'], so m_code needs to be split from name
            pathway = next(temp_pathway for temp_pathway in pathway_list if
                           temp_pathway.map_code == m_code.split(' ',1)[0])  # finds pathway that matches map code
            pathway_name = pathway.name  # stores name of pathway into variable pathway_name
            if target.lower() in pathway_name.lower():  # checks if pathway name matches target name
                    is_there = True
                    break
        return is_there

############ Trims Raw Input File ############
def trim_unannotated(path = None):
    global _input_file
    if path == None:  # if no path is passed to the function
        path = _input_file  # then use the input file path that was previously stored
    if path != _input_file:  # if the input path was different than the _input_file path set
        set_input(path)  # set the new _input_file path and trimmed_file path according to specified path
    try:
        with open(path, 'rU') as gene_list:  # opens specified path of input file to read
            with open(_trimmed_file, 'w+') as output_file:  # creates a new trimmed file write in
                for line in gene_list:
                    CurrGene = line.split()
                    if len(CurrGene) == 2:  # only writes genes that have an ascension number assigned to it
                        output_file.write(line)
            output_file.close()
    except FileNotFoundError:  # if the computer cannot find the path the user specified
        sys.exit("No input file found at " + path)  # exit the program with error message

############ Generate Pathway Map List ############
def gen_pathway(trimmed_file = None, data_file = None):
    global _trimmed_file  # sets local _trimmed_file to global _trimmed_file to use pathway set previously
    global _data_file  # sets local _data_file to global _trimmed_file to use pathway set previously
    if trimmed_file == None:  # if no path to trimmed file is specified,
        trimmed_file = _trimmed_file  # the function will retrieve the trimmed file path that was most recently stored
    if data_file == None:  # if no path to data file is specified,
        if _data_file == "No_File_Specified":  # if _data_file was never set via set_input() or set_data_file(),
            sys.exit("No data file was ever specified.")  # then exit the code with the following msg
        else:
            data_file = _data_file  # the function will retrieve the data file path that was most recently stored
    data_name = os.path.basename(data_file)  # gets the name of the data file without directory and sets to data_name

    if not os.path.isfile(trimmed_file):  # before running rest of code, checks if the trimmed input file exist
        sys.exit("No trimmed input file found at " + trimmed_file)  # if it doesn't exist,exit the program with msg

    pathwayList = []  # creates an empty list to store all the pathways the annotated genes are involved in
    geneList = []  # creates an empty list to store all the annotated genes and their linked pathways
    completed = False  # records whether program has run to completion or not.
    num_already_saved = 0  # records number of gene entries already on saved file
    resume_check_link_path = True  # when resuming, makes sure every gene is added to every pathway map before moving on

    n = 0  # integer to track which gene the code is on

    try:  # trys to open data file
        with open(data_file, "rb") as f:
            if pickle.load(f):  # checks if the data is complete or not
                os.remove(trimmed_file)  # delete the trimmed_file after using it to keep the folder tidy
                sys.exit("A completed data file for " + data_name + " dataset already exist. Please rename or delete " + data_name + " if you wish to work with a rerun dataset.")  # if data is already complete, do not run rest of code
            else:  # if data is not complete, load the data
                num_already_saved = pickle.load(f)  # record where the data left off at
                for _ in range(num_already_saved):  # loads all gene data into geneList
                    geneList.append(pickle.load(f))
                for _ in range(pickle.load(f)):  # loads all pathway data into pathwayList
                    pathwayList.append(pickle.load(f))
    except IOError:
        print("New " + data_name + " will be created.")

    with open(trimmed_file, 'r') as annotated_genes: #opens the trimmed gene list file
        for line in annotated_genes: #iterates through each line of the file
            n = n+1  # for each line, n increases by one to represent it is working on the next code
            print(n)  # print out the number of the gene code is currently working on
            if n > num_already_saved:  # skips all genes that has already been saved
                currGene = line.split()  # isolates gene number from k_code
                gene_number = currGene[0]
                k_code = currGene[1]
                try:
                    geneList.append(Gene(gene_number, k_code))  # creates a gene object from k code and adds it to geneList
                except (ConnectionResetError, TimeoutError) as error:  # if there is a ConnectionResetError or TimeoutError
                    save(geneList, pathwayList, completed)  # save the data to a file
                    sys.exit("Error encountered. Data successfully saved.")  # and exit the program with the following error message
                for path_info in geneList[-1].link_path:  # runs for each pathway the gene was associated with
                    m_code = path_info[:8]  # cuts the map code from the map code + map name of link_path element
                    if any(x.map_code == m_code for x in pathwayList):  # checks to see if pathway was encountered before
                        curr_pathway = next(x for x in pathwayList if x.map_code == m_code)  # finds pathway with that map_code and sets it as curr_pathway
                        if resume_check_link_path:  # check if this is the first pass and need to check which  pathways the gene has already been added to
                            if not any(x == k_code for x in curr_pathway.genes_invol):  # only execute the code directly below if gene hasn't already been added to pathway
                                curr_pathway.add_gene(k_code) # adds gene to that pathway
                        else:  # if this is not the first pass, just run the code directly below
                            curr_pathway.add_gene(k_code)  # t adds gene to that pathway
                    else:  # if this is the first time encountering this pathway
                        try:
                            pathwayList.append(Pathway_MAP(path_info))  # it creates a new pathway object
                        except (ConnectionResetError, TimeoutError) as error:  # if there is a ConnectionResetError or TimeoutError
                            geneList = geneList[:-1]  # deletes the last gene entry from geneList, so when the code restarts it will go through the map_code list of that gene again
                            save(geneList,pathwayList, completed)  # save the data to a file
                            sys.exit("Error encountered. Data successfully saved.")  # and exit the program with the following error maessage
                        pathwayList[-1].add_gene(k_code)  # and adds the gene to the new pathway object
                if n % 100 == 0:  # every 100 genes, run the code directly below
                    save(geneList, pathwayList, completed)  # save the data to a file
            resume_check_link_path = False  # no longer need to check if gene was already added to a pathway after first gene

    completed = True  # if program reaches this step then it is completed
    os.remove(trimmed_file)  # delete the trimmed_file after using it to keep the folder tidy
    save(geneList,pathwayList, completed)  # saves data to file

############ Generate Output File ############
def out_HTML(data_file = None, html_file = None, gene_table = True, pathway_table = True):
    global _html_file  # sets local _html_file to the global _html_file
    global _data_file  # sets local _data_file to the global _data_file
    if data_file == None:  # if no path to data file is passed, then
        data_file = _data_file  # then use the data file path that was previously stored
    else:
        set_data_file(data_file)  # if an argument is passed then, set the _data_file and _html_file paths using argument
    if html_file == None:  # if user does not specify name of html file, then
        html_file = _html_file  # set local html file to the global _html_file

    geneList = _gene_list  # store data of genes in this List
    pathwayList = _pathway_list  # store pathway data in this List

    print("Generating Tables . . . ")  # status update for when program is generating html file

    with open(html_file, "w") as f:
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

        if gene_table:  # if true, includes table for genes
            html_code = html_code + """
    	<div class="scrollingtable">
    		<div>
    			<div>
    				<table>
    					<caption>Genes and Linked Pathway</caption>
    					<thead>
    						<tr>
    			<th><div label="Gene ID"></div></th>
                <th><div label="Gene Name"></div></th>
                <th><div label="K Number"></div></th>
                <th><div label="Definition"></div></th>
    							<th>
    								<!--more versatile way of doing column label; requires 2 identical copies of label-->
    								<div><div>Pathway</div><div>Pathway</div></div>
    							</th>
                  <th class="scrollbarhead"></th> <!--ALWAYS ADD THIS EXTRA CELL AT END OF HEADER ROW-->
    						</tr>
    					</thead>
    					<tbody>
                            """
            for gene in geneList:  # cycles through even gene; one gene each row
                # adds the gene number, k number, and definition to first three columns
                row = "<tr><td>" + gene.gene_num + "</td><td>" + gene.name + "</td><td>" + gene.k_code + "</td><td>" \
                      + gene.definition + "</td><td>"
                first_iteration = True  # used to avoid adding a ", " to the first hyperlink
                for m_code in gene.link_path:  # cycles through each map code for pathways linked to gene
                    m_code = m_code[:8]  # isolates the map code from the map code and name name in link_path item
                    pathway = next(temp_pathway for temp_pathway in pathwayList if
                                   temp_pathway.map_code == m_code)  # finds pathway that matches map code
                    hyper_text = pathway.name + "(" + str(len(
                        pathway.genes_invol)) + ")"  # generates the text for the hyperlink (name + num of associated genes)
                    url = pathway.generate_url(gene.k_code)  # generates the url
                    hyperlink = "<a href=\"" + url + "\">" + hyper_text + "</a>"  # embeds hyperlink to text
                    if not first_iteration:
                        hyperlink = ", " + hyperlink  # adds a ", " between hyperlinks
                    first_iteration = False
                    row = row + hyperlink  # adds hyperlink to row
                row = row + "</td></tr>"  # finishes html code for end of row
                html_code = html_code + row  # add row to html code
            html_code = html_code + """</tbody>
    				</table>
    			</div>
    		</div>
    	</div>"""

        if pathway_table:  # if true, includes the two tables for pathways
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
            pathwayList.sort(key=lambda pathway: len(pathway.genes_invol), reverse=True)
            for pathway in pathwayList:
                hyper_text = pathway.name + "(" + str(
                    len(pathway.genes_invol)) + ")"  # generates the text for the hyperlink (name + num of associated genes)
                url = pathway.generate_url("")  # generates the url
                hyperlink = "<a href=\"" + url + "\">" + hyper_text + "</a>"  # embeds hyperlink to text
                row = "<tr><td>" + hyperlink + "</td></tr>"
                html_code = html_code + row
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
            pathwayList.sort(key=lambda pathway: pathway.name)
            for pathway in pathwayList:
                hyper_text = pathway.name + "(" + str(
                    len(pathway.genes_invol)) + ")"  # generates the text for the hyperlink (name + num of associated genes)
                url = pathway.generate_url("")  # generates the url
                hyperlink = "<a href=\"" + url + "\">" + hyper_text + "</a>"  # embeds hyperlink to text
                row = "<tr><td>" + hyperlink + "</td></tr>"
                html_code = html_code + row
            html_code = html_code + """</tbody>
    				</table>
    			</div>
    		</div>
    	</div>
        </body>
        </html>"""
        f.write(html_code)
    f.close()

    print("Tables Complete")  # status update for when program is finished

def out_txt(output_file = None):
    global _input_file  # sets local _input_file to the global _input_file
    global _html_file  # sets local _html_file to the global _html_file
    global _data_file  # sets local _data_file to the global _data_file

    if output_file == None:  # if user does not specify name of output file, then
        file_name, file_extension = os.path.splitext(_html_file)
        output_file = file_name + "-gene_table" + ".txt"
    output_name, output_extension = os.path.splitext(output_file)  # isolates the output extension from file name
    if output_extension != ".txt":
        output_file = output_name + ".txt"

    geneList = _gene_list  # store data of genes in this List
    pathwayList = _pathway_list  # store pathway data in this List

    with open(output_file, "w") as f:
        f.write("GeneID\tGene_Name\tDefinition\t")
        for gene in geneList:
            line = str(gene.gene_num) + "\t" + str(gene.name) + "\t" + str(gene.definition) + "\t"
            for m_code in gene.link_path:  # cycles through each map code for pathways linked to gene
                m_code = m_code[:8]  # isolates the map code from the map code and name name in link_path item
                pathway = next(temp_pathway for temp_pathway in pathwayList if
                               temp_pathway.map_code == m_code)  # finds pathway that matches map code
                pathway_name = pathway.name  # stores the name of the pathway in a string to manipulate later
                hyper_text = pathway_name.replace(" ", "_") + "(" + str(len(
                    pathway.genes_invol)) + ")"  # generates the text for the hyperlink (name + num of associated genes)
                line = line + hyper_text + " "  # adds pathway name to line to be written along with a space
            line = line + "\n"  # add a new line after finishes storing everything to be written on that line
            f.write(line)

############ User Interface ############
class UI:  # class to wrap all the menu screens that will help user navigate the program
    def menu_data(self):  # first menu that will ask whether to create new data file or use a pre-existing one
        print(textwrap.dedent("""
                                 Would you like to:
                                    1) Create a generate a new data file and table from a new dataset of KO numbers
                                    2) Create a generate a new data file and table from a new dataset of KO numbers (Batch)
                                    3) Generate a new table from an existing data file
                                    4) Generate a new table from an existing data file (Batch)
                                """))  # displays options for
        choice = input("Input a digit for your choice: ")  # ask user for input as a single digit
        if choice == '1':  # if user chose '1'
            self.menu_data_new()  # then initiate menu branch for creating data file and html table from new dataset
        elif choice == '2':  # if user chose '2'
            self.menu_batch_list(data = False)  # then initiate menu branch for processing multiple input files at once
        elif choice == '3':  # if user chose '3'
            self.menu_data_existing()  # then initiate menu branch that just generates html table from pre-existing data file
        elif choice == '4':  # if user chose '4'
            self.menu_batch_list(data = True)  # then initiate menu branch for processing multiple data files at once
        else:  # if input was not 1 or 2, then ask again
            print("Not a valid choice")
            self.menu_data()


    def menu_data_new(self, input_file = None, input_list = None):  # menu that prompts for input file to create a new data file from dataset of KO numbers
        print(textwrap.dedent("""
                                 Please enter either the relative or absolute path to the input file below"""))
        if input_file == None:
            input_file = input()  # accepts user specified path to input file and stores as string in input_file

        set_input(input_file)  # uses user specified path to set paths and file names the code will use
        if os.path.isfile(_input_file):  # checks if the file specified by the user exists
            print(textwrap.dedent("""
                                     Trimming input file"""))  # status update
            trim_unannotated()  # trim the data file of genes that are not associated with a KO number
            print(textwrap.dedent("""
                                     Trimming input file (complete)
                                     
                                     Extracting data from KEGG"""))  #status update
            gen_pathway()  # accesses KEGG API to extract information on genes and pathways associated with input KO#s
            load_data()  # loads data from newly generated data file to populate _gene_list and _pathway_list

            print(textwrap.dedent("""Extracting data from KEGG (complete)"""))  # status update
            self.menu_filters(input_list = input_list)
        else:  # if file could not be found, then re-prompt for input file location
            print(_input_file + " was not found")
            self.menu_data_new(input_list = input_list)

    def menu_data_existing(self, data_file = None, input_list = None):  # does not extract any info from KEGG, but uses pre-generated data file
        print(textwrap.dedent("""
                                 Please enter either the relative or absolute path to the data file below
                                 """))
        if data_file == None:
            data_file = input()  # accepts user specified path to data file and stores as string in data_file
        file_name, extension = os.path.splitext(data_file)  # stores file name and extension respectively
        if extension == '':  # if no extension is found...
            data_file = data_file + ".dat"  # add the .dat extension onto the path
            extension = ".dat"  # sets extension to ".dat" to pass the next if statement
        if extension.lower() == ".dat":  # checks if the file has the correct extension
            try:  # try to open data file
                with open(data_file, 'rb') as f:
                    if pickle.load(f):  # loads first value from data file which is a boolean indicating if file is complete
                        set_data_file(data_file)  # if file is complete, then set paths using the path the user input
                        load_data(data_file)  # loads lists in data file into global variables to be used
                        self.menu_filters(input_list = input_list)
                    else:  # if the file is incomplete, return to first menu
                        print("Data file is incomplete")
                        self.menu_data()  # returns user to first menu where they may complete the data file
            except FileNotFoundError:  # if file could not be found, then re-prompt for data file location
                print(data_file + " was not found")
                self.menu_data_existing(input_list = input_list)
        else:  # if file does not have the correct format ".dat"
            print(textwrap.dedent("""
                                         Please choose a dat file as the data file."""))  # inform user of error
            self.menu_data_existing(input_list = input_list)  # reprompt the question

    def menu_batch_list(self, data):  # menu that prompts for list of input files to create a new data file from dataset of KO numbers
        # data = True for using data files; data = False when using input files
        file_list = []  # list of input/data file paths
        if data:
            print(textwrap.dedent("""
                                     Please enter either the relative or absolute path to a text file containing a list of the data files.
                                     Example text file with list of data files:
                                     C_elegans-1.dat
                                     C_elegans-2.dat
                                     D_pulex-1.dat

                                     Enter the word 'all' to use every dat file in the current directory as data files.
                                     Ensure all data files being run are complete
                                     """))
        else:
            print(textwrap.dedent("""
                                     Please enter either the relative or absolute path to a text file containing a list of the input files.
                                     Example text file with list of input files:
                                     C_elegans-1.txt
                                     C_elegans-2.txt
                                     D_pulex-1.txt
                                     
                                     Enter the word 'all' to use every txt file in the current directory as input files.
                                     Ensure all files being run are input txt files and that they do NOT have a complete data file already.
                                     """))
        choice = input()  # accepts user specified path to text file with list of input/data files; may also accept "all"
        if choice.lower() == 'all':  # if user wants to use all .txt/.dat files in current directory
            file_list = retrieve_batch_files(data)  # stores all txt/dat files in current directory in file_list
            self.menu_batch_run(data, file_list)  # run a batch run using the specified file list
        else:  # if user specifies a path to a text file containing a list of input/data files instead
            try:  # try to open the text file
                with open(choice, 'r') as f:  # opens the text file in order to extract input file paths
                    for line in f:  # cycle through each line of the list of input files
                        file_list.append(line.rstrip())  # add line to file_list after stripping the \n away
                self.menu_batch_run(data, file_list)  # run a batch run using the specified file list
            except FileNotFoundError:
                print(choice + " was not found")
                self.menu_batch_list(data)

    def menu_batch_run(self, data, file_list):
        user_inputs = []  # list to store user inputs for batch runs
        user_inputs.extend(self.menu_filters(batch_ask = True))  # stores user's filter choices in user_input list
        if user_inputs[0]  == '1':  # if user decided to add a filter to the data
            user_inputs.append('2')  # add choice of '2) No' when prompted to add another filter
        user_inputs.append("")  # adds choice of no custom name
        user_inputs.append(self.menu_table_type(batch_ask = True, custom_name = None))  # adds choice of table type

        for file in file_list:
            temp_inputs = user_inputs[:]
            if data:
                self.menu_data_existing(data_file = file, input_list = temp_inputs)
            else:
                self.menu_data_new(input_file = file, input_list = temp_inputs)

    def menu_filters(self, batch_ask = False, input_list = None):  # ask whether user would like to filter data
        filter_present = len(_gene_list) != _total_genes  # determines whether a filter is present
        if not filter_present:
            print(textwrap.dedent("""
                                     Would you like to filter the data?
                                        1) Yes
                                        2) No
                                     """))
        else:
            print(textwrap.dedent("""
                                     Would you like to apply another filter?
                                        1) Yes
                                        2) No
                                        3) Remove all filters
                                     """))
        if input_list:  # if this is a batch run with a list of user inputs available...
            choice = input_list.pop(0)  # use the choice from the user input list instead of prompting a new input
        else:  # if no input list is provided
            choice = input()  # prompt the user to choose whether or not to filter data
        if choice == '1':
            if batch_ask:  # if running this function only to record user input for batch run
                filter_choices = [choice]  # stores all of user choice of using a filter
                filter_choices.extend(self.menu_filters_type(batch_ask))  # stores user choice of filter type and search
                if filter_choices[-1] == '0':  # if user changed mind about using a filter
                    filter_choices = ['2']  # change first choice to '2) No' to save time
                return filter_choices  # return user choices for filter
            self.menu_filters_type(input_list = input_list)
        elif choice == '2':
            if batch_ask:  # if running this function only to record user input for batch run
                return choice  # return user input to not filter data
            self.menu_table_name(input_list = input_list)
        elif choice == '3':
            load_data()  # program will repopulate gene and pathway lists to reverse effects of filters
            self.menu_filters(input_list = input_list)
        else:
            print("Not a valid choice")
            self.menu_filters(batch_ask, input_list = input_list)

    def menu_filters_type(self, batch_ask = False, input_list = None):  # asks user which type of filter they would like to apply to the data
        global _gene_list #  no longer refers to a local global _gene_list
        print(textwrap.dedent("""
                                 What would you like to filter by?
                                    1) Name
                                    2) Definition
                                    3) Pathway
                                    4) Go back to previous menu
                                 """))

        if input_list:  # if this is a batch run with a list of user inputs available...
            filter_type = input_list.pop(0)  # use the choice from the user input list instead of prompting a new input
        else:
            filter_type = input()
        valid_choice = True  # keeps track of whether a valid input was made

        if filter_type != '4':  # if user wants to go back to previous menu w/o filtering then skip below code

            print(textwrap.dedent("""
                                     Enter text you would like to search for: 
                                     """))
            if input_list:  # if this is a batch run with a list of user inputs available...
                search = input_list.pop(0)  # use the choice from the user input list instead of prompting a new input
            else:
                search = input()

            if batch_ask:  # if running this function only to record user input for batch run
                return [filter_type, search]  # return user input for filter type and search term

            prev_gene_list = _gene_list[:]

            if filter_type == '1':
                print("Number of genes before filter: " + str(len(_gene_list)))  # displays number of genes before filter
                # only keep genes that return True for check_name
                _gene_list[:] = [gene for gene in _gene_list if gene.check_name(search)]
                print("Number of genes after filter: " + str(len(_gene_list)))  # displays num of genes after filter
            elif filter_type == '2':
                print("Number of genes before filter: " + str(len(_gene_list)))  # displays number of genes before filter
                # only keep genes that return True for check_definition
                _gene_list[:] = [gene for gene in _gene_list if gene.check_definition(search)]
                print("Number of genes after filter: " + str(len(_gene_list)))  # displays num of genes after filter
            elif filter_type == '3':
                print("Number of genes before filter: " + str(len(_gene_list)))  # displays number of genes before filter
                # only keep genes that return True for check_pathway
                _gene_list[:] = [gene for gene in _gene_list if gene.check_pathway(search, _pathway_list)]
                print("Number of genes after filter: " + str(len(_gene_list)))  # displays num of genes after filter
            else:
                valid_choice = False  # if input was not 1-3, then the input was not valid, so valid_choice is False

            if len(_gene_list) == 0:
                print("No entries found with that search. Reverting filter.")
                _gene_list = prev_gene_list[:]
                print("Number of genes after previous filter reverted: " + str(len(_gene_list)))

        if valid_choice:  # test to see whether a valid choice was made
            if batch_ask:  # if chooses not to apply filter afterall during questions for batch run...
                return '0'  # return 0 to alert program that no filter will be used (batch run can only ask once currently)
            else:
                self.menu_filters()  # if a valid choice was made then go back to the menu.filters menu
        else:
            print("Not a valid choice")  # if a valid choice was not made, then print
            self.menu_filters_type(input_list = input_list)  # and return to same menu for allow re-input

    def menu_table_name(self, input_list = None):  # this menu prompts user to choose a name for the html file
        custom_name = None  # might store user defined html file name
        basename, ext = get_file_name(_html_file)  # cuts path from file name
        name_no_ext, ext = os.path.splitext(basename)  # stores file name w/o ext in name_no_ext and the extension in ext
        print(textwrap.dedent("""
                                 Enter the output file name or press ENTER to use the default name [""" + name_no_ext + """]:
                                 """))
        if input_list:  # if this is a batch run with a list of user inputs available...
            choice = input_list.pop(0)  # use the choice from the user input list instead of prompting a new input
        else:
            choice = input()  # prompts user to set name of output file
        if choice.strip() != "":
            custom_name = os.path.join(os.path.dirname(_html_file), choice + ext)  # creates absolute path from specified name
        self.menu_table_type(custom_name = custom_name, input_list = input_list)

    def menu_table_type(self, custom_name, batch_ask = False, input_list = None):  # this menu prompts user to choose between making a table for genes or pathways

        print(textwrap.dedent("""
                                 Would you like to:
                                    1) Generate a table of genes (HTML)
                                    2) Generate a table of genes and pathways (HTML)
                                    3) Generate a table of pathways (HTML)
                                    4) Generate a table of genes without links to pathway maps (txt, tab-delimited, small size)
                                 """))
        if input_list:  # if this is a batch run with a list of user inputs available...
            choice = input_list.pop(0)  # use the choice from the user input list instead of prompting a new input
        else:
            choice = input("Input a digit for your choice: ")

        if batch_ask:  # if running this function only to record user input for batch run
            return choice  # return user's choice of table type

        print("\nCreating table\n")  # status update so user knows that script is processing

        # if pathway tables needs to be generated and some genes have been filtered
        if (choice == '2' or choice == '3') and len(_gene_list) != _total_genes:
            # trims _pathway_list to remove pathways where with no associated genes (were removed in filter)
            _pathway_list[:] = [pathway for pathway in _pathway_list if pathway.check_genes(_gene_list)]

        if choice == '1':
            out_HTML(html_file = custom_name, pathway_table = False)  # generates only the genes table
        elif choice == '2':
            out_HTML(html_file = custom_name)  # generates both genes and pathways tables, and run next menu to set html file name
        elif choice == '3':
            out_HTML(html_file = custom_name, gene_table = False)  # generates only the pathways table
        elif choice == '4':
            out_txt(output_file = custom_name)
        else:
            print("Not a valid choice")
            self.menu_table_type()

############ if script is run then do this ############
if __name__ == "__main__":
    ui = UI()
    ui.menu_data()
