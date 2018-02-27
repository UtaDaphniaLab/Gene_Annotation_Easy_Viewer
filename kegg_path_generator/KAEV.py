############ Imports ############
import contextlib  # the closing context manager will ensure connection to url is closed after with statement
from urllib.request import urlopen  # import the tools we need to open url
import bisect  # module needed to insert an item into a sorted list
import pickle  # needed to save and load data
import sys  # needed to exit out of program when error is encountered
import re  # required to split files without removing delimiter

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
    with open("PathGeneData.dat", "wb") as f:  # create a file to save/pickle data
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

class Pathway_MAP:  # class for pathway map objects
    def __init__(self, ipathway_info):  # accepts the map code of the pathway
        pathway_info = ipathway_info.split(None, 1)  # separtes the map code and the map name by spiting at first space
        self.map_code = pathway_info[0]  # sotres map code in map_code
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
        info_url = decode_url(link)  # access KEGG for information
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

############ Trims Raw Input File ############
def trim_unannotated():
    input_file = input("Path to input file: ")  # ask user for absolute or relative path to input file
    try:
        with open(input_file, 'rU') as gene_list:  # opens specified path of input file to read
            in_file_sep = separate_file(input_file)  # will separate the file name from rest of the path
            if in_file_sep[0]:   # simply adds "trimmed_" in front of old name for new file
                output_file_name = str(in_file_sep[0]) + "trimemed_" + str(in_file_sep[1])  # use for absolute path
            else:
                output_file_name = "trimmed_" + str(in_file_sep[1])  # use for relative path
            with open(output_file_name, 'w+') as output_file:  # creates a new trimmed file write in
                CurrGene = []  # makes a list to store the current gene the program is examining

                for line in gene_list:
                    CurrGene = line.split()
                    if len(CurrGene) == 2:  # only writes genes that have an ascension number assigned to it
                        output_file.write(line)
            output_file.close()
    except FileNotFoundError:  # if the computer cannot find the path the user specified
        sys.exit("No file found at \"" + input_file + "\"")  # exit the program with error message
    return output_file_name

############ Generate Pathway Map List ############
def gen_pathway(trimmed_input_file):
    pathwayList = []  # creates an empty list to store all the pathways the annotated genes are involved in
    geneList = []  # creates an empty list to store all the annotated genes and their linked pathways
    completed = False  # records whether program has run to completion or not.
    num_already_saved = 0  # records number of gene entries already on saved file
    resume_check_link_path = True  # when resuming, makes sure every gene is added to every pathway map before moving on

    n = 0  # integer to track which gene the code is on

    try:  # trys to open PathGeneData.dat
        with open("PathGeneData.dat", "rb") as f:
            if pickle.load(f):  # checks if the data is complete or not
                sys.exit("A completed data file already exist. Please rename PathGeneData.dat if you wish to work with a different dataset.")  # if data is already complete, do not run rest of code
            else:  # if data is not complete, load the data
                num_already_saved = pickle.load(f)  # record where the data left off at
                for _ in range(num_already_saved):  # loads all gene data into geneList
                    geneList.append(pickle.load(f))
                for _ in range(pickle.load(f)):  # loads all pathway data into pathwayList
                    pathwayList.append(pickle.load(f))
        f.closed
    except IOError:
        print("New PathGeneData.dat will be created.")

    with open(trimmed_input_file, 'r') as annotated_genes: #opens the trimed gene list file
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
    save(geneList,pathwayList, completed)  # saves data to file

############ Generate HTML Output File ############
def out_HTML():
    geneList = []  # store data of genes in this List
    pathwayList = []  # store pathway data in this List

    with open("PathGeneData.dat", "rb") as f:  # open data file that was generated previously
        if not pickle.load(f):  # reads the first line of saved data which tells whether it is complete or not
            sys.exit("Data is not complete")  # exits program if data is not complete

        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in genes
            geneList.append(pickle.load(f))  # loads data into geneList
        for _ in range(pickle.load(f)):  # reads line that tells program how many data entries there are in pathways
            pathwayList.append(pickle.load(f))  # loads data into pathwayList

    with open("Annotated Gene Pathways.html", "w") as f:
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
        <body link="#003366">
    	<div class="scrollingtable">
    		<div>
    			<div>
    				<table>
    					<caption>Genes and Linked Pathway</caption>
    					<thead>
    						<tr>
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
            row = "<tr><td>" + gene.name + "</td><td>" + gene.k_code + "</td><td>" + gene.definition + "</td><td>"
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
        html_code = html_code + """					</tbody>
    				</table>
    			</div>
    		</div>
    	</div>
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
        html_code = html_code + """					</tbody>
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
        html_code = html_code + """					</tbody>
    				</table>
    			</div>
    		</div>
    	</div>
        </body>
        </html>"""
        f.write(html_code)
    f.close()

############ if script is run then do this ############
if __name__ == "__main__":
    trimmed_file_name = trim_unannotated()
    gen_pathway(trimmed_file_name)
    out_HTML()
