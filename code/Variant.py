from flask import Flask,render_template,flash,request,url_for,redirect,jsonify,Response,send_from_directory
import csv
import MySQLdb
import MySQLdb.cursors
import json
import re
from decimal import Decimal
import vcf
import os
from werkzeug import secure_filename
import time


UPLOAD_FOLDER='/home/bioinfo/Heart_gene_database/HeartInter/Uploads/'
ALLOWED_EXTENSIONS=set(['txt','vcf'])
app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

delimiters = '\t',';',"\r\n",',',' '
regexPattern = '|'.join(map(re.escape,delimiters))

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]),re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

"""file reader"""
def allowed_file(filename):
    return '.' in filename and \
            filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS


"""get allele frequency from database"""
def get_freq(con,table,column,poplen,chr,pos,ref,alt):
    """Connect to MySQL database"""
    sql_command="SELECT %s FROM allele_frequency.%s WHERE (chr =\"%s\") AND (pos=\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (column,table,chr,pos,ref,alt)
    con.execute(sql_command)
    freq_result = {}
    result = con.fetchall()
    if len(result) > 0:
        freq_result = dict(result[0])
    else:
        from itertools import repeat
        value = []
        value.extend(repeat(".",poplen))
        freq_result = dict(zip(list(column.split(",")),value))
    return freq_result


@app.route("/")
def index():
    return render_template('index.html')

@app.route("/search",methods=['GET','POST'])
def search():
    if request.method == "POST":
#connecting to the database
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Expression_profiles",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

#processing user input gene list
        Search_text=multiple_replace(request.form['search_list'],{"\"":"","\'":"","`":"","%":""})
        Search_text=Search_text.upper()
        Search_list =list(filter(None,re.split(regexPattern,Search_text)))

#the annotation table user choose. Currently, one allow one species
        Annotation_table = request.form["data_Annotation"] +"_"+request.form.get('ref')+"_RNA"
        Annotation_title = ["User_input","ensembl_gene_id","gene_name_ensembl","gene_name_NCBI","aliases","chr","gene_start","gene_end","strand","gene_description","family_ID"]

        Annotation_column = ",".join(list(map(lambda orig_string:Annotation_table+"."+orig_string,Annotation_title[1:])))

        #tissue list user select
        Human_exp = request.form.getlist("Human_exp")
        Mouse_exp = request.form.getlist("Mouse_exp")
        Zebrafish_exp = request.form.getlist("Zebrafish_exp")

        if len(Human_exp)!=0:
            if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                Annotation_title.append("Human_ID")
                Annotation_title.extend(Human_exp)
        if len(Mouse_exp)!=0:
            if "Mouse_GRCm38_RNA" != Annotation_table:
                Annotation_title.append("Mouse_ID")
                Annotation_title.extend(Mouse_exp)
        if len(Zebrafish_exp)!=0:
            if "Zebrafish_GRCz10_RNA" != Annotation_table:
                Annotation_title.append("Zebrafish_ID")
                Annotation_title.extend(Zebrafish_exp)

        species = set()
        species.add(Annotation_table)


#title processing
        exp_table = request.form["data_Annotation"] + "_exp"
        exp_column = ",".join(list(map(lambda orig_string:Annotation_table+"."+orig_string+"_FPKM"+","+Annotation_table+"."+orig_string +"_Rank",eval(exp_table))))
        Annotation_title.extend(exp_column.replace(Annotation_table +".","").split(","))
        Annotation_column = ",".join([Annotation_column,exp_column])

        FinalResult = {}
        Dupl = set()
        NoData = set()
        FinalResultList = []
        if request.form.get('SearchType') == "Symbol":
            Ensembl_Set = set()
            NCBI_Set = set()
            aliase_Set = set()
            Ensembl_Gene = 'SELECT %s FROM `%s` WHERE `gene_name_ensembl` IN (%s)'
            in_p = ','.join(list(map(lambda x:'\'' + x + '\'', Search_list)))
            Ensembl_Gene = Ensembl_Gene % (Annotation_column,Annotation_table,in_p)
            con.execute(Ensembl_Gene)
            Result = con.fetchall()
            for Ensembl_Res in Result:
                if Ensembl_Res["gene_name_ensembl"].upper() in FinalResult:
                    Ensembl_Res["User_input"] = Ensembl_Res["gene_name_ensembl"].upper()
                    FinalResult[Ensembl_Res["gene_name_ensembl"].upper()].append(Ensembl_Res)
                    Dupl.add(Ensembl_Res["gene_name_ensembl"].upper())


                else:
                    FinalResult[Ensembl_Res["gene_name_ensembl"].upper()] = []
                    Ensembl_Res["User_input"] = Ensembl_Res["gene_name_ensembl"].upper()
                    FinalResult[Ensembl_Res["gene_name_ensembl"].upper()].append(Ensembl_Res)
                    Ensembl_Set.add(Ensembl_Res["gene_name_ensembl"].upper())

            if len(set(Search_list)-Ensembl_Set)!=0:
                in_p = ', '.join(list(map(lambda x: '\'' + x +'\'',set(Search_list)-Ensembl_Set)))
                NCBI_Gene = 'SELECT %s FROM `%s` WHERE `gene_name_NCBI` IN (%s)'
                NCBI_Gene = NCBI_Gene % (Annotation_column,Annotation_table,in_p)
                con.execute(NCBI_Gene)
                Result = con.fetchall()
                for NCBI_Res in Result:
                    if NCBI_Res["gene_name_NCBI"].upper() in FinalResult:
                        NCBI_Res["User_input"] = NCBI_Res["gene_name_NCBI"].upper()
                        FinalResult[NCBI_Res["gene_name_NCBI"].upper()].append(NCBI_Res)
                        Dupl.add(NCBI_Res["gene_name_NCBI"].upper())
                    else:
                        FinalResult[NCBI_Res["gene_name_NCBI"].upper()] = []
                        NCBI_Res["User_input"] = NCBI_Res["gene_name_NCBI"].upper()
                        FinalResult[NCBI_Res["gene_name_NCBI"].upper()].append(NCBI_Res)
                        NCBI_Set.add(NCBI_Res["gene_name_NCBI"].upper())
            aliase = set(Search_list) - Ensembl_Set - NCBI_Set
            if len(set(aliase))!=0:
                for symbol in aliase:
                    aliase_Gene = 'SELECT %s FROM `%s` WHERE FIND_IN_SET("%s",`aliases`)' % (Annotation_column,Annotation_table,symbol)
                    con.execute(aliase_Gene)
                    Result = con.fetchall()
                    if len(Result) != 0:
                        Result = [dict(Raliase, User_input=symbol) for Raliase in Result]
                        aliase_Set.add(symbol)
                        FinalResult[symbol]=[]
                        FinalResult[symbol].extend(list(Result))
                        if len(Result)>1:
                            Dupl.add(symbol)
            NoData = aliase - aliase_Set


        elif request.form.get('SearchType') == "Ensembl_id":
            Ensembl_Set = set()
            Ensembl_Gene = 'SELECT %s FROM `%s` WHERE `ensembl_gene_id` IN (%s)'
            in_p = ','.join(list(map(lambda x:'\'' + x + '\'',Search_list)))
            Ensembl_Gene = Ensembl_Gene % (Annotation_column,Annotation_table,in_p)
            con.execute(Ensembl_Gene)
            Result = con.fetchall()
            for Ensembl_Res in Result:
                if Ensembl_Res["ensembl_gene_id"].upper() in FinalResult:
                    Ensembl_Res["User_input"] =Ensembl_Res["ensembl_gene_id"].upper()
                    FinalResult[Ensembl_Res["ensembl_gene_id"].upper()].append(Ensembl_Res)
                    Dupl.add(Ensembl_Res["ensembl_gene_id"].upper())
                else:
                    FinalResult[Ensembl_Res["ensembl_gene_id"].upper()] = []
                    Ensembl_Res["User_input"] =Ensembl_Res["ensembl_gene_id"].upper()
                    FinalResult[Ensembl_Res["ensembl_gene_id"].upper()].append(Ensembl_Res)
                    Ensembl_Set.add(Ensembl_Res["ensembl_gene_id"].upper())
            NoData = set(Search_list)-Ensembl_Set



        for genename,values in FinalResult.items():
            if genename not in Dupl:
                if len(Human_exp)!=0:
                    if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                        species.add("Human_GRCh37_RNA")
                        sql_command = 'SELECT  GROUP_CONCAT(CONCAT_WS(":",ensembl_gene_id,%s)) AS %s , %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                        sql_command = sql_command % \
                        (",".join(list(map(lambda orig_string:"\""+orig_string+"\","+orig_string+"_FPKM"+","+orig_string+"_Rank",Human_exp))),"Human_ID",",".join(list(map(lambda orig_string:"AVG("+orig_string+"_FPKM) AS "+orig_string,Human_exp))),"Human_GRCh37_RNA",values[0]['family_ID'])
                        con.execute(sql_command)
                        Result_exp = con.fetchall()
                        values[0].update(Result_exp[0])
                        FinalResultList.append(values[0])


                if len(Mouse_exp)!=0:
                    if "Mouse_GRCm38_RNA" != Annotation_table:
                        species.add("Mouse_GRCm38_RNA")
                if len(Zebrafish_exp)!=0:
                    if "Zebrafish_GRCz10_RNA" != Annotation_table:
                        species.add("Zebrafish_column")
                        sql_command = 'SELECT GROUP_CONCAT(CONCAT_WS(":",ensembl_gene_id,%s)) AS %s , %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                        sql_command = sql_command % \
                        (",".join(list(map(lambda orig_string:"\""+orig_string+"\","+orig_string+"_FPKM"+","+orig_string+"_Rank",Zebrafish_exp))),"Zebrafish_ID",",".join(list(map(lambda orig_string:"AVG("+orig_string+"_FPKM) AS "+ orig_string,Zebrafish_exp))),"Zebrafish_GRCz10_RNA",values[0]['family_ID'])
                        con.execute(sql_command)
                        Result_exp = con.fetchall()
                        values[0].update(Result_exp[0])
                        FinalResultList.append(values[0])

        if len(NoData)!=0:
            for no in NoData:
                value=[no]
                from itertools import repeat
                value.extend(list(repeat(".",len(Annotation_title)-1)))
                NoData_Res = dict(zip(Annotation_title,value))
                FinalResultList.append(NoData_Res)
        db.close()

        if len(Dupl)!=0:
            return render_template("geneselect.html",results=FinalResult,keys=Dupl,Final=FinalResultList,Annotation_title=Annotation_title,Annotation_table=Annotation_table,Human_exp=Human_exp,Mouse_exp=Mouse_exp,Zebrafish_exp=Zebrafish_exp)

        return render_template("results.html",results=FinalResultList,keys=Annotation_title)
    return render_template('new_search.html')

@app.route('/geneselect', methods=['GET','POST'])
def geneselect():
    if request.method == 'POST':

#connecting to the database
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Expression_profiles",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

        dupl = request.form["duplname"][1:-1].replace("\'","").replace(" ","").split(",")
        FinalString= request.form["Final"].replace("\'","\"")
        Final=json.loads(FinalString)
        Annotation_title = request.form["Annotation_title"]

        Annotation_title =json.loads(request.form["Annotation_title"].replace("\'","\""))
        Annotation_table = request.form["Annotation_table"]
        Human_exp = json.loads(request.form["Human_exp"].replace("\'","\""))
        Zebrafish_exp = json.loads(request.form["Zebrafish_exp"].replace("\'","\""))
        print(Annotation_table)

        for duplgene in dupl:
            dup_res = json.loads(request.form[duplgene].replace("\'","\""))
            if len(Human_exp)!=0:
                if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                    sql_command = 'SELECT GROUP_CONCAT(CONCAT_WS(":",ensembl_gene_id,%s)) AS %s , %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                    sql_command = sql_command % (",".join(list(map(lambda orig_string:"\""+orig_string+"\","+orig_string+"_FPKM"+","+orig_string+"_Rank",Human_exp))),"Human_ID",",".join(list(map(lambda orig_string:"AVG("+orig_string+"_FPKM) AS "+orig_string,Human_exp))),"Human_GRCh37_RNA",dup_res['family_ID'])
                    con.execute(sql_command)
                    Result_exp = con.fetchall()
                    dup_res.update(Result_exp[0])
            if len(Zebrafish_exp)!=0:
                if "Zebrafish_GRCz10_RNA" != Annotation_table:
                    sql_command = 'SELECT GROUP_CONCAT(CONCAT_WS(":",ensembl_gene_id,%s)) AS %s , %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                    sql_command = sql_command % (",".join(list(map(lambda orig_string:"\""+orig_string+"\","+orig_string+"_FPKM"+","+orig_string+"_Rank",Zebrafish_exp))),"Zebrafish_ID",",".join(list(map(lambda orig_string:"AVG("+orig_string+"_FPKM) AS "+ orig_string,Zebrafish_exp))),"Zebrafish_GRCz10_RNA",dup_res['family_ID'])
                    con.execute(sql_command)
                    Result_exp = con.fetchall()
                    dup_res.update(Result_exp[0])
            Final.append(dup_res)
            print(Final)
            db.close()
        return render_template("results.html",results=Final,keys=Annotation_title)
    return render_template("geneselect.html")


@app.route('/vcf', methods=['GET','POST'])
def upload_file():

    if request.method == 'POST':
        """Create a uuid"""
        import uuid
        User_id = str(uuid.uuid1())

        from itertools import repeat

        """Connect to MySQL database"""
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

        """user input information"""
        Output_format = request.form["format"]
        Genomes_population = request.form.getlist("Genomes")#1000 Genomes population list
        JPN_population = request.form.getlist("JPN")
        ESP_population = request.form.getlist("ESP")
        tissue = request.form.getlist("tissue")
        REVEL_threshold = float(request.form["REVEL_threshold"])


        """table title"""
        Final_result = [] #存最後的結果
        Final_result_title = ["chr","pos","ref","alt"]

        """add user select tissue to gene annotation list"""
        Gene_annotation =["chr","pos","ref","alt","gene_name_ensembl","gene_description"]

        tissue_column = ",".join(list(map(lambda orig_string:orig_string +"_FPKM"+","+orig_string + "_Rank",tissue)))
        Gene_annotation.extend(tissue_column.split(","))

       #add different population to allele freq table

        Population_allele_freq = ["chr","pos","ref","alt"]
        Predict=["chr","pos","ref","alt","Index"]
        if len(Genomes_population)!=0:
            Genomes_population_column = ",".join(list(map(lambda orig_string:orig_string + "_Ref_" + Output_format +","+orig_string + "_Alt_"+ Output_format,Genomes_population)))
            Final_result_title.extend(Genomes_population_column.split(","))
            Population_allele_freq.extend(Genomes_population_column.split(","))
        if "1KJPN" in JPN_population:
            JPN_1_column = "1KJPN_Ref_" + Output_format + "," + "1KJPN_Alt_" +Output_format
            Final_result_title.extend(JPN_1_column.split(","))
            Population_allele_freq.extend(JPN_1_column.split(","))
        if "2KJPN" in JPN_population:
            JPN_2_column = "2KJPN_Ref_" + Output_format + "," + "2KJPN_Alt_" +Output_format
            Final_result_title.extend(JPN_2_column.split(","))
            Population_allele_freq.extend(JPN_2_column.split(","))
        if len(ESP_population) !=0:
            ESP_population_column = ",".join(list(map(lambda orig_string:orig_string + "_Ref_" +Output_format +","+orig_string + "_Alt_" +Output_format,ESP_population)))
            Final_result_title.extend(ESP_population_column.split(","))
            Population_allele_freq.extend(ESP_population_column.split(","))
        if len(request.form.getlist("TWB"))!=0:
            TWB_column = "TWB_Ref_" + Output_format + "," + "TWB_Alt_" +Output_format
            Final_result_title.extend(TWB_column.split(","))
            Population_allele_freq.extend(TWB_column.split(","))
        if len(request.form.getlist("REVEL")) != 0:
            Final_result_title.extend(["aaref","aaalt","REVEL"])
            Predict.extend(["aaref","aaalt","REVEL"])
        Final_result_title.extend(["gene_name_ensembl","gene_description"])

        Final_result_title.extend(tissue_column.split(","))
        Predict.extend(["gerp++","Func","ExonicFunc","AAChange"])

        Result_line = 0


        #if user input file
        file = request.files['file']
        #如果使用者上傳檔案，則使用此檔案的資訊
        if file.filename == '':
            variant_text = multiple_replace(request.form['variants_list'],{"\"":"","\'":"","`":"","%":""})
            variant_text = variant_text.upper()
            variant_list =list(filter(None,re.split(regexPattern,variant_text)))
            vartext_file = open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".avinput"),"w")

            for var in variant_list:
                vars = re.compile("([\dXYM]+):(\d+)([ATCG]+)>([ATCG]+)").split(var)
                Final_result.append({})
                if len(vars)<6:
                    Final_result[Result_line].update({"chr":vars[0],"pos":".","ref":".","alt":"."})
                    vartext_file.write(var[0]+"\t.\t.\t.\n")
                else:
                    Final_result[Result_line].update({"chr":vars[1],"pos":vars[2],"ref":vars[3],"alt":vars[4]})
                    vartext_file.write("chr"+vars[1]+"\t"+vars[2]+"\t"+str(int(vars[2])+len(vars[3])-1)+"\t"+vars[3]+"\t"+vars[4]+"\n")

                if len(Genomes_population)!=0:
                    sql_command="SELECT %s FROM allele_frequency.1000Genomes_5pop_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (Genomes_population_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    Genomes_result = con.fetchall()
                    if len(Genomes_result)>0:
                        Final_result[Result_line].update(Genomes_result[0])
                    else:
                        value = []
                        value.extend(list(repeat(".",len(Genomes_population)*2)))
                        Genomes_result = dict(zip(list(Genomes_population_column.split(",")),value))
                        Final_result[Result_line].update(Genomes_result)



                if "1KJPN" in JPN_population:
                    sql_command="SELECT %s FROM allele_frequency.1KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (JPN_1_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    JPN1_result = con.fetchall()
                    if len(JPN1_result)>0:
                        Final_result[Result_line].update(JPN1_result[0])
                    else:

                        value=[]
                        value.extend(list(repeat(".",2)))
                        JPN1_result =dict(zip(list(JPN_1_column.split(",")),value))
                        Final_result[Result_line].update(JPN1_result)

                if "2KJPN" in JPN_population:
                    sql_command="SELECT %s FROM allele_frequency.2KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (JPN_2_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    JPN2_result = con.fetchall()
                    if len(JPN2_result)>0:
                        Final_result[Result_line].update(JPN2_result[0])
                    else:
                        value=[]
                        value.extend(repeat(".",2))
                        JPN2_result =dict(zip(list(JPN_2_column.split(",")),value))
                        Final_result[Result_line].update(JPN2_result)

                if len(ESP_population) !=0:
                    sql_command="SELECT %s FROM allele_frequency.ESP_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (ESP_population_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    ESP_result = con.fetchall()
                    if len(ESP_result)>0:
                        Final_result[Result_line].update(ESP_result[0])
                    else:
                        value=[]
                        value.extend(repeat(".",len(ESP_population)*2))
                        ESP_result =dict(zip(list(ESP_population_column.split(",")),value))
                        Final_result[Result_line].update(ESP_result)
                if len(request.form.getlist("TWB"))!=0:
                    sql_command="SELECT %s FROM allele_frequency.TWB_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt=\"%s\")" % (TWB_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    TWB_result = con.fetchall()
                    if len(TWB_result)>0:
                        Final_result[Result_line].update(TWB_result[0])
                    else:
                        value=[]
                        value.extend(repeat(".",2))
                        TWB_result =dict(zip(list(TWB_column.split(",")),value))
                        Final_result[Result_line].update(TWB_result)




                if len(tissue)!=0:
                    sql_command = "SELECT `gene_name_ensembl`, `gene_description`,H_heart_muscle_Rank,%s FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (tissue_column,vars[1],vars[2])
                if len(tissue)==0:
                    sql_command = "SELECT `gene_name_ensembl`,`gene_description`,`H_heart_muscle_Rank` FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (vars[1],vars[2])

                con.execute(sql_command)
                RNA_result = con.fetchall()
                if len(RNA_result) > 0:
                    Final_result[Result_line].update(RNA_result[0])
                else:
                    from itertools import repeat
                    Final_result[Result_line].update({"gene_name_ensembl":".","gene_description":"."})
                    value=[]
                    value.extend(repeat(".",len(tissue)*2))
                    RNA_result =dict(zip(list(tissue_column.split(",")),value))
                    Final_result[Result_line].update(RNA_result)
                Final_result[Result_line].update(get_freq(con,"REVEL","aaref,aaalt,REVEL",3,vars[1],vars[2],vars[3],vars[4]))
                Result_line +=1
            vartext_file.close()

        elif file and allowed_file(file.filename):
            filename = User_id + secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
            vcf_reader = vcf.Reader(open(os.path.join(app.config['UPLOAD_FOLDER'],filename)))

            for record in vcf_reader:
                for idx in range(0,len(record.ALT)):
                    Final_result.append({})
                    Final_result[Result_line].update({"chr":record.CHROM,"pos":record.POS,"ref":record.REF,"alt":record.ALT[idx]})
                    if len(Genomes_population)!=0:
                        sql_command="SELECT %s FROM allele_frequency.1000Genomes_5pop_%s WHERE (chr=\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (Genomes_population_column,Output_format,record.CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        Genomes_result = con.fetchall()
                        if len(Genomes_result)>0:
                            Final_result[Result_line].update(Genomes_result[0])
                        else:
                            value = []
                            value.extend(repeat(".",len(Genomes_population)*2))
                            Genomes_result = dict(zip(list(Genomes_population_column.split(",")),value))
                            Final_result[Result_line].update(Genomes_result)

                    if "1KJPN" in JPN_population:
                        sql_command="SELECT %s FROM allele_frequency.1KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (JPN_1_column,Output_format,record.CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        JPN1_result = con.fetchall()
                        if len(JPN1_result)>0:
                            Final_result[Result_line].update(JPN1_result[0])
                        else:
                            value=[]
                            value.extend(repeat(".",2))
                            JPN1_result =dict(zip(list(JPN_1_column.split(",")),value))
                            Final_result[Result_line].update(JPN1_result)
                    if "2KJPN" in JPN_population:
                        sql_command="SELECT %s FROM allele_frequency.2KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (JPN_2_column,Output_format,record.CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        JPN2_result = con.fetchall()
                        if len(JPN2_result)>0:
                            Final_result[Result_line].update(JPN2_result[0])
                        else:
                            value=[]
                            value.extend(repeat(".",2))
                            JPN2_result =dict(zip(list(JPN_2_column.split(",")),value))
                            Final_result[Result_line].update(JPN2_result)
                    if len(ESP_population) !=0:
                        sql_command="SELECT %s FROM allele_frequency.ESP_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (ESP_population_column,Output_format,record.CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        ESP_result = con.fetchall()
                        if len(ESP_result)>0:
                            Final_result[Result_line].update(ESP_result[0])
                        else:
                            value=[]
                            value.extend(repeat(".",len(ESP_population)*2))
                            ESP_result =dict(zip(list(ESP_population_column.split(",")),value))
                            Final_result[Result_line].update(ESP_result)
                    if len(request.form.getlist("TWB"))!=0:
                        sql_command="SELECT %s FROM allele_frequency.TWB_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (TWB_column,Output_format,record.CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        TWB_result = con.fetchall()
                        if len(TWB_result)>0:
                            Final_result[Result_line].update(TWB_result[0])
                        else:
                            value=[]
                            value.extend(repeat(".",2))
                            TWB_result =dict(zip(list(TWB_column.split(",")),value))
                            Final_result[Result_line].update(TWB_result)

                    #RNA expression
                    if len(tissue)!=0:
                        sql_command = "SELECT `gene_name_ensembl`,`gene_description`,H_heart_muscle_Rank,%s FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (tissue_column,record.CHROM,record.POS)
                    if len(tissue)==0:
                        sql_command = "SELECT `gene_name_ensembl`,`gene_description`,H_heart_muscle_Rank FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (record.CHROM,record.POS)
                    con.execute(sql_command)
                    RNA_result = con.fetchall()
                    if len(RNA_result) > 0:
                        Final_result[Result_line].update(RNA_result[0])
                    else:
                        Final_result[Result_line].update({"gene_name_ensembl":".","gene_description":"."})
                        value=[]
                        value.extend(repeat(".",len(tissue)*2))
                        RNA_result =dict(zip(list(tissue_column.split(",")),value))
                        Final_result[Result_line].update(RNA_result)
                    Final_result[Result_line].update(get_freq(con,"REVEL","aaref,aaalt,REVEL",3,record.CHROM,record.POS,record.REF,record.ALT[idx]))
                    Result_line +=1

       #ANNOVAR

        import subprocess
        import sys


        if file.filename == '':
            cmd='/home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/table_annovar.pl %s /home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/humandb/ -buildver hg19 -remove -protocol ensGene,gerp++gt2 -operation g,f -nastring . --outfile %s' % (os.path.join(app.config['UPLOAD_FOLDER'],User_id +".avinput"),os.path.join(app.config['UPLOAD_FOLDER'],User_id))

        else:
            cmd='/home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/table_annovar.pl %s /home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/humandb/ -buildver hg19 -remove -protocol ensGene,gerp++gt2 -operation g,f -nastring . -vcfinput --outfile %s' % (os.path.join(app.config['UPLOAD_FOLDER'],filename),os.path.join(app.config['UPLOAD_FOLDER'],User_id))

        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0: sys.exit(retcode)
        Current_line = 0
        Final_result_title.extend(["Func.refGene","ExonicFunc.refGene","AAChange.refGene","gerp++gt2"])
        with open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".hg19_multianno.txt"),"r") as annovar:
            for line in annovar:
                Current_line+=1
                line = line.strip()
                if Current_line !=1:
                    lines = line.split("\t")
                    Final_result[Current_line-2].update({"Func":lines[5],"ExonicFunc":lines[8],"AAChange":lines[9],"gerp++":lines[10]})

        db.close()
        """Add index"""
        for resaddID in Final_result:
            if resaddID['Func'] == "intergenic" or resaddID['Func'] == "intronic" or resaddID['ExonicFunc']=="synonymous SNV":
                resaddID['Index'] = 0
            if resaddID['REVEL']!= ".":
                if resaddID['ExonicFunc'] == "nonsynonymous SNV" or float(resaddID['REVEL'])<REVEL_threshold:
                    resaddID['Index'] = 1
                if resaddID['ExonicFunc'] == "nonsynonymous SNV":
                    if float(resaddID['REVEL'])>REVEL_threshold or resaddID['Func'] == "splicing":
                        resaddID['Index'] = 2
                        if resaddID['H_heart_muscle_Rank']!=".":
                            if eval(resaddID['H_heart_muscle_Rank'])<0.25:
                                resaddID['Index'] = 3
                        if resaddID['gerp++']!=".":
                                resaddID['Index'] = 3
            elif resaddID['REVEL']== ".":
                if resaddID['ExonicFunc'] == "nonsynonymous SNV":
                    resaddID['Index'] = 1
                if resaddID['Func'] == "splicing":
                    resaddID['Index'] = 2
                    if resaddID['H_heart_muscle_Rank']!=".":
                        if eval(resaddID['H_heart_muscle_Rank'])<0.25:
                            resaddID['Index'] = 3
                    if resaddID['gerp++']!=".":
                            resaddID['Index'] = 3
            if "Index" not in resaddID:
                resaddID['Index'] = 0


        """remove user files"""
        import glob
        for rmfile in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], User_id+"*")):
            os.remove(rmfile)

        return render_template("Teresult.html",results = Final_result,Gene_annotation=Gene_annotation,Population_allele_freq=Population_allele_freq,Predict=Predict)
    return render_template('vcf.html')


@app.route("/Tutorial")
def Tutorial():
    return render_template('Tutorial.html')


if __name__=='__main__':
    app.run(host='0.0.0.0',port=5000,debug=True,threaded=True)
