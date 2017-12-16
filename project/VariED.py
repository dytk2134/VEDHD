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
from itertools import repeat
from collections import namedtuple
import gzip

#path
#apache

UPLOAD_FOLDER='/var/www/html/varied/upload/'
InterVar_Path='/var/www/html/varied/annovar/Intervar.py'
CADD_path = '/home/bioinfo/Variedtools/CADD_v1.3/bin/score.sh'
ANNOVAR_path = '/var/www/html/varied/annovar/table_annovar.pl'
humandb_path = '/var/www/html/varied/annovar/humandb/'

#test
"""
UPLOAD_FOLDER='/home/bioinfo/Heart_gene_database/HeartInter/Uploads/'
InterVar_Path='/var/www/html/varied/annovar/Intervar.py'
CADD_path = '/home/bioinfo/Variedtools/CADD_v1.3/bin/score.sh'
ANNOVAR_path = '/var/www/html/varied/annovar/table_annovar.pl'
humandb_path = '/var/www/html/varied/annovar/humandb/'
"""
ALLOWED_EXTENSIONS=set(['txt','vcf'])
app = Flask(__name__)
app.secret_key = "^awed@1qh)#1ozd0+2dx*d117l3cr!@rfnr238jducwmpt0cd_"
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

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/Expression_profiles",methods=['GET','POST'])
def Expression_profiles():
    if request.method == "POST":
#connecting to the database"""
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Expression_profiles",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

#processing user input gene list"""
        Search_text=multiple_replace(request.form['search_list'],{"\"":"","\'":"","`":"","%":""})
        Search_text=Search_text.upper()
        Search_list =list(filter(None,re.split(regexPattern,Search_text)))

#the annotation table user choose. Currently, only allow one species"""
        Annotation_table = request.form["data_Annotation"] +"_"+request.form.get('ref')+"_RNA"
        Annotation_title = ["User_input","ensembl_gene_id","gene_name_ensembl","gene_name_NCBI","aliases","chr","gene_start","gene_end","strand","gene_description"]

        Annotation_column = ",".join(list(map(lambda orig_string:Annotation_table+"."+orig_string,Annotation_title[1:])))
        Annotation_column = Annotation_column + ",family_ID"


        #tissue list user select
        Human_exp = request.form.getlist("Human_exp")
        Mouse_exp = request.form.getlist("Mouse_exp")
        Zebrafish_exp = request.form.getlist("Zebrafish_exp")
        Human_orthologs = ["gene_name_ensembl","ensembl_gene_id"]
        Mouse_orthologs = ["gene_name_ensembl","ensembl_gene_id"]
        Zebrafish_orthologs = ["gene_name_ensembl","ensembl_gene_id"]

#title processing
        exp_table = request.form["data_Annotation"] + "_exp"
        exp_column = ",".join(list(map(lambda orig_string:Annotation_table+"."+orig_string+"_FPKM"+","+Annotation_table+"."+orig_string +"_Rank",eval(exp_table))))
        Annotation_title.extend(exp_column.replace(Annotation_table +".","").split(","))
        Annotation_column = ",".join([Annotation_column,exp_column])
        if Annotation_column.endswith(','):
            Annotation_column = Annotation_column[:-1]
        if len(Human_exp)!=0:
            if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                Annotation_title.append("Human_orthologs")

                Annotation_title.extend([s + "_AVGExp" for s in Human_exp])
                Orthologs = ",".join(list(map(lambda orig_string:orig_string+"_FPKM"+","+orig_string +"_Rank",Human_exp)))
                Human_orthologs.extend(Orthologs.split(","))

        if len(Mouse_exp)!=0:
            if "Mouse_GRCm38_RNA" != Annotation_table:
                Annotation_title.append("Mouse_orthologs")
                Annotation_title.extend([s + "_AVGExp" for s in Mouse_exp])
                Orthologs = ",".join(list(map(lambda orig_string:orig_string+"_FPKM"+","+orig_string +"_Rank",Mouse_exp)))
                Mouse_orthologs.extend(Orthologs.split(","))

        if len(Zebrafish_exp)!=0:
            if "Zebrafish_GRCz10_RNA" != Annotation_table:
                Annotation_title.append("Zebrafish_orthologs")
                Annotation_title.extend([s + "_AVGExp" for s in Zebrafish_exp])
                Orthologs = ",".join(list(map(lambda orig_string:orig_string+"_FPKM"+","+orig_string +"_Rank",Zebrafish_exp)))
                Zebrafish_orthologs.extend(Orthologs.split(","))


        species = set()
        species.add(Annotation_table)


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
                if in_p.endswith(","):
                    in_p = in_p[:-1]
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
                if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                    if len(Human_exp)!=0:
                        if values[0]['family_ID']!=None:
                            sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                            sql_command = sql_command % \
                        (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+orig_string+"_AVGExp",Human_exp))),"Human_GRCh37_RNA",values[0]['family_ID'])
                            con.execute(sql_command)
                            Result_exp = con.fetchall()
                            if len(Result_exp)!=0:
                                values[0].update(Result_exp[0])
                                sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Human_GRCh37_RNA WHERE `family_ID` = \"%s\"'
                                sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Human_exp))),values[0]['family_ID'])
                                con.execute(sql_command)
                                values[0]["Human_orthologs"] = list(con.fetchall())
                        if "Human_orthologs" not in values[0]:
                            value=[]
                            value.extend(list(repeat("None",len(Human_exp))))
                            Result_exp =dict(zip([s + "_AVGExp" for s in Human_exp],value))
                            values[0].update(Result_exp)
                            values[0]["Human_orthologs"] = "None"

                if "Mouse_GRCm38_RNA" != Annotation_table:
                    if len(Mouse_exp)!=0:
                        if values[0]['family_ID']!=None:
                            sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                            sql_command = sql_command % \
                        (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+orig_string+"_AVGExp",Mouse_exp))),"Mouse_GRCm38_RNA",values[0]['family_ID'])
                            con.execute(sql_command)
                            Result_exp = con.fetchall()
                            if len(Result_exp)!=0:
                                values[0].update(Result_exp[0])
                                sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Mouse_GRCm38_RNA WHERE `family_ID` = \"%s\"'
                                sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Mouse_exp))),values[0]['family_ID'])
                                con.execute(sql_command)
                                values[0]["Mouse_orthologs"] = list(con.fetchall())
                        if "Mouse_orthologs" not in values[0]:
                            value=[]
                            value.extend(list(repeat("None",len(Mouse_exp))))
                            Result_exp =dict(zip([s + "_AVGExp" for s in Mouse_exp],value))
                            values[0].update(Result_exp)
                            values[0]["Mouse_orthologs"] = "None"

                if "Zebrafish_GRCz10_RNA" != Annotation_table:
                    if len(Zebrafish_exp)!=0:
                        if values[0]['family_ID']!=None:
                            sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                            sql_command = sql_command % \
                            (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+ orig_string+"_AVGExp",Zebrafish_exp))),"Zebrafish_GRCz10_RNA",values[0]['family_ID'])
                            con.execute(sql_command)
                            Result_exp = con.fetchall()
                            if len(Result_exp)!=0:
                                values[0].update(Result_exp[0])
                                sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Zebrafish_GRCz10_RNA WHERE `family_ID` = \"%s\"'
                                sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Zebrafish_exp))),values[0]['family_ID'])
                                con.execute(sql_command)
                                values[0]["Zebrafish_orthologs"] = list(con.fetchall())
                        if "Zebrafish_orthologs" not in values[0]:
                            value=[]
                            value.extend(list(repeat("None",len(Zebrafish_exp))))
                            Result_exp =dict(zip([s + "_AVGExp" for s in Zebrafish_exp],value))
                            values[0].update(Result_exp)
                            values[0]["Zebrafish_orthologs"] = "None"

                FinalResultList.append(values[0])

        if len(NoData)!=0:
            for no in NoData:
                value=[no]
                value.extend(list(repeat("None",len(Annotation_title)-1)))
                NoData_Res = dict(zip(Annotation_title,value))
                FinalResultList.append(NoData_Res)

        db.close()

        if len(Dupl)!=0:
            return render_template("geneselect.html",results=FinalResult,keys=Dupl,Final=FinalResultList,Annotation_title=Annotation_title,Annotation_table=Annotation_table,Human_exp=Human_exp,Mouse_exp=Mouse_exp,Zebrafish_exp=Zebrafish_exp,Human_orthologs=Human_orthologs,Mouse_orthologs=Mouse_orthologs,Zebrafish_orthologs=Zebrafish_orthologs)


        return render_template("Expression_profiles_result.html",results=FinalResultList,keys=Annotation_title,Human_orthologs=Human_orthologs,Mouse_orthologs=Mouse_orthologs,Zebrafish_orthologs=Zebrafish_orthologs,Annotation_table=Annotation_table)
    return render_template('Expression_profiles.html')

@app.route('/geneselect', methods=['GET','POST'])
def geneselect():
    if request.method == 'POST':

#connecting to the database
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Expression_profiles",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

        dupl = request.form["duplname"][1:-1].replace("\'","").replace(" ","").split(",")
        FinalString= request.form["Final"].replace("\'","\"").replace("None","null")
        Final=json.loads(FinalString)
        Annotation_title = request.form["Annotation_title"]

        Annotation_title =json.loads(request.form["Annotation_title"].replace("\'","\""))
        Annotation_table = request.form["Annotation_table"]
        Human_exp = json.loads(request.form["Human_exp"].replace("\'","\""))
        Mouse_exp= json.loads(request.form["Mouse_exp"].replace("\'","\""))
        Zebrafish_exp = json.loads(request.form["Zebrafish_exp"].replace("\'","\""))

        Human_orthologs= json.loads(request.form["Human_orthologs"].replace("\'","\""))
        Mouse_orthologs= json.loads(request.form["Mouse_orthologs"].replace("\'","\""))
        Zebrafish_orthologs= json.loads(request.form["Zebrafish_orthologs"].replace("\'","\""))


        for duplgene in dupl:
            dup_res = json.loads(request.form[duplgene].replace("\'","\"").replace("None","null"))
            if "Human_GRCh37_RNA" != Annotation_table and "Human_GRCh38_RNA"!= Annotation_table:
                if len(Human_exp)!=0:
                    if dup_res['family_ID']!=None:
                        sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                        sql_command = sql_command % (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+orig_string+"_AVGExp",Human_exp))),"Human_GRCh37_RNA",dup_res['family_ID'])
                        con.execute(sql_command)
                        Result_exp = con.fetchall()
                        if len(Result_exp)!=0:
                            dup_res.update(Result_exp[0])
                            sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Human_GRCh37_RNA WHERE `family_ID` = \"%s\"'
                            sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Human_exp))),dup_res['family_ID'])
                            con.execute(sql_command)
                            dup_res["Human_orthologs"] = list(con.fetchall())
                    if "Human_orthologs" not in dup_res:
                        value=[]
                        value.extend(list(repeat("None",len(Human_exp))))
                        Result_exp =dict(zip([s + "_AVGExp" for s in Human_exp],value))
                        dup_res.update(Result_exp)
                        dup_res["Human_orthologs"] = "None"

            if "Zebrafish_GRCz10_RNA" != Annotation_table:
                if len(Zebrafish_exp)!=0:
                    if dup_res['family_ID']!=None:
                        sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                        sql_command = sql_command % (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+ orig_string+"_AVGExp",Zebrafish_exp))),"Zebrafish_GRCz10_RNA",dup_res['family_ID'])
                        print(sql_command)
                        con.execute(sql_command)
                        Result_exp = con.fetchall()
                        if len(Result_exp)!=0:
                            dup_res.update(Result_exp[0])
                            sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Zebrafish_GRCz10_RNA WHERE `family_ID` = \"%s\"'
                            sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Zebrafish_exp))),dup_res['family_ID'])
                            con.execute(sql_command)
                            dup_res["Zebrafish_orthologs"] = list(con.fetchall())
                    if "Zebrafish_orthologs" not in dup_res:
                        value=[]
                        value.extend(list(repeat("None",len(Zebrafish_exp))))
                        Result_exp =dict(zip([s + "_AVGExp" for s in Zebrafish_exp],value))
                        dup_res.update(Result_exp)
                        dup_res["Zebrafish_orthologs"] = "None"

            if "Mouse_GRCm38_RNA" != Annotation_table:
                if len(Mouse_exp)!=0:
                    if dup_res['family_ID']!=None:
                        sql_command = 'SELECT %s FROM %s WHERE `family_ID` = \"%s\" GROUP BY family_ID'
                        sql_command = sql_command % (",".join(list(map(lambda orig_string:"ROUND(AVG("+orig_string+"_FPKM),2) AS "+ orig_string+"_AVGExp",Mouse_exp))),"Mouse_GRCm38_RNA",dup_res['family_ID'])
                        con.execute(sql_command)
                        Result_exp = con.fetchall()
                        if len(Result_exp)!=0:
                            dup_res.update(Result_exp[0])
                            sql_command = 'SELECT gene_name_ensembl,ensembl_gene_id,%s FROM Mouse_GRCm38_RNA WHERE `family_ID` = \"%s\"'
                            sql_command = sql_command % (",".join(list(map(lambda orig_string:orig_string+"_FPKM,"+orig_string+"_Rank",Mouse_exp))),dup_res['family_ID'])
                            con.execute(sql_command)
                            dup_res["Mouse_orthologs"] = list(con.fetchall())
                    if "Mouse_orthologs" not in dup_res:
                        value=[]
                        value.extend(list(repeat("None",len(Mouse_exp))))
                        Result_exp =dict(zip([s + "_AVGExp" for s in Mouse_exp],value))
                        dup_res.update(Result_exp)
                        dup_res["Mouse_orthologs"] = "None"

            Final.append(dup_res)

        return render_template("Expression_profiles_result.html",results=Final,keys=Annotation_title,Human_orthologs=Human_orthologs,Mouse_orthologs=Mouse_orthologs,Zebrafish_orthologs=Zebrafish_orthologs)
    return render_template("geneselect.html")

User_id =""
@app.route('/Variants_search', methods=['GET','POST'])
def Variants_search():
    if request.method == 'POST':

        """Create a uuid"""
        import uuid
        User_id =str(uuid.uuid1())
        """Connect to MySQL database"""
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()

        """user input information"""
        Output_format = request.form["format"]
        Genomes_population = request.form.getlist("Genomes")#1000 Genomes population list
        JPN_population = request.form.getlist("JPN")
        ESP_population = request.form.getlist("ESP")
        TaiwanBiobank = request.form.getlist("TWB")
        tissue = request.form.getlist("tissue")
        REVEL_threshold = float(request.form["REVEL_threshold"])
        detected_level = 0.5

        """table title"""
        Final_result = [] #存最後的結果
        Final_result_title = ["chr","pos","ref","alt","gene_name_ensembl","gene_description"]

        """add user select tissue to gene annotation list"""
        Gene_annotation =["chr","pos","ref","alt","gene_name_ensembl","gene_description"]
        tissue_column_search = ",".join(list(map(lambda orig_string:orig_string +"_FPKM"+","+orig_string + "_Rank",tissue)))
        tissue_column_out = ""
        exp_out_list = request.form.getlist("exp_out")
        exp_out_list.append("detected")
        for tis in tissue:
            for exp_out in exp_out_list:
                tissue_column_out += "," + tis +"_"+exp_out
        if len(tissue)!=0:

            Gene_annotation.extend(tissue_column_out.split(",")[1:])
            Final_result_title.extend(tissue_column_out.split(",")[1:])
        Gene_annotation.append("External")

       #add different population to allele freq table

        Population_allele_freq = ["chr","pos","ref","alt"]
        Predict=["chr","pos","ref","alt","Index"]
        Clinical_interpretation = ["chr","pos","ref","alt"]
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
            if "TWB_GWG" in TaiwanBiobank:
                TWB_GWG_column = "TWB_GWG_Ref_" + Output_format + "," + "TWB_GWG_Alt_" +Output_format
                Final_result_title.extend(TWB_GWG_column.split(","))
                Population_allele_freq.extend(TWB_GWG_column.split(","))
            if "TWB_NGS" in TaiwanBiobank:
                TWB_NGS_column = "TWB_NGS_Ref_" + Output_format + "," + "TWB_NGS_Alt_" +Output_format
                Final_result_title.extend(TWB_NGS_column.split(","))
                Population_allele_freq.extend(TWB_NGS_column.split(","))
        if len(request.form.getlist("REVEL")) != 0:
            Final_result_title.extend(["REVEL"])
            Predict.extend(["REVEL"])
        if len(request.form.getlist("InterVar")) !=0:
            Final_result_title.append("InterVar")
            Clinical_interpretation.append("InterVar")
        if len(request.form.getlist("ClinVar"))!=0:
            Final_result_title.extend(["AlleleID","ClinicalSignificance","PhenotypeList"])
            Clinical_interpretation.extend(["AlleleID","ClinicalSignificance","PhenotypeList"])
        if len(request.form.getlist("CADD"))!=0:
            Final_result_title.extend(["CADD_RawScore","CADD_PHRED"])
            Predict.extend(["CADD_RawScore","CADD_PHRED"])
        target_tissue = request.form["target_tissue"]

        ANNOVAR_column = ['Func.ensGene','ExonicFunc.ensGene','AAChange.ensGene','gerp++gt2']
        Final_result_title.extend(["Index","gerp++gt2","Func.ensGene","ExonicFunc.ensGene","AAChange.ensGene"])
        Predict.extend(["gerp++gt2","Func.ensGene","ExonicFunc.ensGene","AAChange.ensGene"])
        if "gnomad_exome" in request.form.getlist("anno"):
            ANNOVAR_column.extend(request.form.getlist("gnomad_exome"))
            Final_result_title.extend(request.form.getlist("gnomad_exome"))
            Population_allele_freq.extend(request.form.getlist("gnomad_exome"))
        if "gnomad_genome" in request.form.getlist("anno"):
            ANNOVAR_column.extend(request.form.getlist("gnomad_genome"))
            Final_result_title.extend(request.form.getlist("gnomad_genome"))
            Population_allele_freq.extend(request.form.getlist("gnomad_genome"))
            print(request.form.getlist("gnomad_genome"))
        if "dbnsfp30a" in request.form.getlist("anno"):
            ANNOVAR_column.extend(request.form.getlist("Polyphen2"))
            Final_result_title.extend(request.form.getlist("Polyphen2"))
            Predict.extend(request.form.getlist("Polyphen2"))


        Result_line = 0
        line_dict = {}#for CADD result


        #if user input file
        file = request.files['file']
        #如果使用者上傳檔案，則使用此檔案的資訊
        if file.filename == '':
            variant_text = multiple_replace(request.form['variants_list'],{"\"":"","\'":"","`":"","%":""})
            variant_text = variant_text.upper()
            variant_list =list(filter(None,re.split(regexPattern,variant_text)))
            vartext_file = open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".avinput"),"w")
            if len(request.form.getlist("CADD"))!=0:
                CADD_file =gzip.open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".vcf.gz"),"wb")

            if len(variant_list)==0:
                flash('Please input at least one variant.')
                removeuserfile(User_id)
                vartext_file.close()
                if len(request.form.getlist("CADD"))!=0:
                    CADD_file.close()

                return redirect(url_for('Variants_search'))

            for var in variant_list:

                vars = re.compile("([\dXYM]+):(\d+)([ATCG-]+)>([ATCG-]+)").split(var)
                Final_result.append({})

                if len(vars)<6 or len(vars)>6 or (len(vars)==6 and vars[5]!=""):
                    Final_result[Result_line].update({"chr":var,"pos":".","ref":".","alt":".","REVEL":".","External":"."})
                    vartext_file.write(var+"\t.\t.\t.\t.\n")
                    value = []
                    value.extend(list(repeat(".",len(Final_result_title[4:]))))
                    Notfound = dict(zip(Final_result_title[4:],value))
                    Final_result[Result_line].update(Notfound)
                    Result_line+=1
                    continue
                else:
                    Final_result[Result_line].update({"chr":vars[1],"pos":vars[2],"ref":vars[3],"alt":vars[4]})
                    vartext_file.write("chr"+vars[1]+"\t"+vars[2]+"\t"+str(int(vars[2])+len(vars[3])-1)+"\t"+vars[3]+"\t"+vars[4]+"\n")
                    if len(request.form.getlist("CADD"))!=0:
                        CADD_line = vars[1]+"\t"+vars[2]+"\t"+".\t"+vars[3]+"\t"+vars[4]+"\n"
                        CADD_file.write(bytes(CADD_line, 'UTF-8'))
                        line_dict_key = (vars[1],vars[2],vars[3],vars[4])
                        if line_dict_key not in line_dict:
                            line_dict[line_dict_key] = []
                        line_dict[line_dict_key].append(Result_line)


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
                if "TWB_GWG" in TaiwanBiobank:
                    sql_command="SELECT %s FROM allele_frequency.TWB_GWG_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt=\"%s\")" % (TWB_GWG_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    TWB_result = con.fetchall()
                    if len(TWB_result)>0:
                        Final_result[Result_line].update(TWB_result[0])
                    else:
                        value=[]
                        value.extend(repeat(".",2))
                        TWB_result =dict(zip(list(TWB_GWG_column.split(",")),value))
                        Final_result[Result_line].update(TWB_result)
                if "TWB_NGS" in TaiwanBiobank:
                    sql_command="SELECT %s FROM allele_frequency.TWB_NGS_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt=\"%s\")" % (TWB_NGS_column,Output_format,vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    TWB_result = con.fetchall()
                    if len(TWB_result)>0:
                        Final_result[Result_line].update(TWB_result[0])
                    else:
                        value=[]
                        value.extend(repeat(".",2))
                        TWB_result =dict(zip(list(TWB_NGS_column.split(",")),value))
                        Final_result[Result_line].update(TWB_result)



                if len(tissue)!=0:
                    sql_command = "SELECT `ensembl_gene_id`,`gene_name_ensembl`, `gene_description`,`%s`,%s FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (target_tissue+"_FPKM",tissue_column_search,vars[1],vars[2])
                if len(tissue)==0:
                    sql_command = "SELECT `ensembl_gene_id`,`gene_name_ensembl`,`gene_description`,`%s` FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (target_tissue+"_FPKM",vars[1],vars[2])

                con.execute(sql_command)
                RNA_result = con.fetchall()
                if len(RNA_result) > 0:
                    Final_result[Result_line].update(RNA_result[0])
                    for RNA_k,RNA_v in RNA_result[0].items():
                        if "FPKM" in RNA_k:
                            if RNA_v != "." and RNA_v !=None:
                                if float(str(RNA_v)) >= detected_level:
                                    Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"Yes"})
                                else:
                                    Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"No"})
                            else:
                                Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"."})

                    URL = "."
                    if RNA_result[0]["ensembl_gene_id"]!=".":
                        URL ="<a href=\"http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="+RNA_result[0]["ensembl_gene_id"]+"\">Ensembl</a>"
                    if RNA_result[0]["gene_name_ensembl"]!=".":
                        URL =URL+"," + "<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term="+RNA_result[0]["gene_name_ensembl"]+"\">NCBI</a>"
                        URL =URL+"," + "<a href=\"https://ghr.nlm.nih.gov/search?query="+RNA_result[0]["gene_name_ensembl"]+"\">NIH - GHI</a>"
                    Final_result[Result_line].update({"External":URL})
                else:
                    Final_result[Result_line].update({"gene_name_ensembl":".","gene_description":".","ensembl_gene_id":".","External":"."})
                    value=[]
                    value.extend(repeat(".",len(tissue_column_out.split(","))-1))
                    RNA_result =dict(zip(list(tissue_column_out.split(","))[1:],value))
                    Final_result[Result_line].update(RNA_result)
                sql_command="SELECT %s FROM allele_frequency.%s WHERE (chr =\"%s\") AND (pos=\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % ("REVEL","REVEL",vars[1],vars[2],vars[3],vars[4])
                con.execute(sql_command)
                result = con.fetchall()
                if len(result) > 0:
                    Final_result[Result_line].update(result[0])
                else:
                    Final_result[Result_line].update({"REVEL":"."})
                if len(request.form.getlist("ClinVar"))!=0:
                    sql_command = "SELECT `AlleleID`,`ClinicalSignificance`,`PhenotypeList` FROM allele_frequency.ClinVar WHERE (Chromosome =\"%s\") AND (Start=\"%s\") AND (ReferenceAllele = \"%s\") AND (AlternateAllele = \"%s\")" % (vars[1],vars[2],vars[3],vars[4])
                    con.execute(sql_command)
                    result = con.fetchall()
                    if len(result) > 0:
                        Final_result[Result_line].update(result[0])
                    else:
                        value = []
                        value.extend(repeat(".",3))
                        ClinVar_result = dict(zip(["AlleleID","ClinicalSignificance","PhenotypeList"],value))
                        Final_result[Result_line].update(ClinVar_result)
                if len(request.form.getlist("CADD"))!=0:
                    Final_result[Result_line].update({"CADD_RawScore":".","CADD_PHRED":"." })



                Result_line +=1
            vartext_file.close()
            if len(request.form.getlist("CADD"))!=0:
                CADD_file.close()


        elif file and allowed_file(file.filename):
            filename = User_id + secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
            if len(request.form.getlist("CADD"))!=0:
                CADD_file =gzip.open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".vcf.gz"),"wb")

            try:
                vcf_reader = vcf.Reader(open(os.path.join(app.config['UPLOAD_FOLDER'],filename)))
                for record in vcf_reader:
                    for idx in range(0,len(record.ALT)):
                        Final_result.append({})
                        CHROM = str(record.CHROM).upper().replace("CHR","")
                        Final_result[Result_line].update({"chr":CHROM,"pos":record.POS,"ref":record.REF,"alt":record.ALT[idx]})
                        if len(request.form.getlist("CADD"))!=0:
                            CADD_line = "%s\t%s\t.\t%s\t%s\n" % (CHROM,record.POS,record.REF,record.ALT[idx])
                            CADD_file.write(bytes(CADD_line, 'UTF-8'))
                            line_dict_key = (CHROM,str(record.POS),str(record.REF),str(record.ALT[idx]))

                            if line_dict_key not in line_dict:
                                line_dict[line_dict_key] = []
                            line_dict[line_dict_key].append(Result_line)


                        if len(Genomes_population)!=0:
                            sql_command="SELECT %s FROM allele_frequency.1000Genomes_5pop_%s WHERE (chr=\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (Genomes_population_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
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
                            sql_command="SELECT %s FROM allele_frequency.1KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (JPN_1_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
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
                            sql_command="SELECT %s FROM allele_frequency.2KJPN_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (JPN_2_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
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
                            sql_command="SELECT %s FROM allele_frequency.ESP_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (ESP_population_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
                            con.execute(sql_command)
                            ESP_result = con.fetchall()
                            if len(ESP_result)>0:
                                Final_result[Result_line].update(ESP_result[0])
                            else:
                                value=[]
                                value.extend(repeat(".",len(ESP_population)*2))
                                ESP_result =dict(zip(list(ESP_population_column.split(",")),value))
                                Final_result[Result_line].update(ESP_result)
                        if "TWB_GWG" in TaiwanBiobank:
                            sql_command="SELECT %s FROM allele_frequency.TWB_GWG_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (TWB_GWG_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
                            con.execute(sql_command)
                            TWB_result = con.fetchall()
                            if len(TWB_result)>0:
                                Final_result[Result_line].update(TWB_result[0])
                            else:
                                value=[]
                                value.extend(repeat(".",2))
                                TWB_result =dict(zip(list(TWB_GWG_column.split(",")),value))
                                Final_result[Result_line].update(TWB_result)
                        if "TWB_NGS" in TaiwanBiobank:
                            sql_command="SELECT %s FROM allele_frequency.TWB_NGS_%s WHERE (chr =\"%s\") AND (pos =\"%s\") AND (ref = \"%s\") AND (alt= \"%s\")" % (TWB_NGS_column,Output_format,CHROM,record.POS,record.REF,record.ALT[idx])
                            con.execute(sql_command)
                            TWB_result = con.fetchall()
                            if len(TWB_result)>0:
                                Final_result[Result_line].update(TWB_result[0])
                            else:
                                value=[]
                                value.extend(repeat(".",2))
                                TWB_result =dict(zip(list(TWB_NGS_column.split(",")),value))
                                Final_result[Result_line].update(TWB_result)





                        #RNA expression
                        if len(tissue)!=0:
                            sql_command = "SELECT `ensembl_gene_id`,`gene_name_ensembl`,`gene_description`,%s,%s FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (target_tissue+"_FPKM",tissue_column_search,CHROM,record.POS)
                        if len(tissue)==0:
                            sql_command = "SELECT `ensembl_gene_id`,`gene_name_ensembl`,`gene_description`,%s FROM Expression_profiles.Human_GRCh37_RNA WHERE chr = \"%s\" AND %s BETWEEN Human_GRCh37_RNA.gene_start AND Human_GRCh37_RNA.gene_end" % (target_tissue+"_FPKM",CHROM,record.POS)
                        con.execute(sql_command)
                        RNA_result = con.fetchall()
                        if len(RNA_result) > 0:
                            Final_result[Result_line].update(RNA_result[0])
                            for RNA_k,RNA_v in RNA_result[0].items():
                                if "FPKM" in RNA_k:
                                    if RNA_v != "." and RNA_v !=None:
                                        if float(str(RNA_v)) >= detected_level:
                                            Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"Yes"})
                                        else:
                                            Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"No"})
                                    else:
                                        Final_result[Result_line].update({RNA_k.replace("_FPKM","_detected"):"."})
                            URL="."
                            if RNA_result[0]["ensembl_gene_id"]!=".":
                                URL ="<a href=\"http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g="+RNA_result[0]["ensembl_gene_id"]+"\">Ensembl</a>"
                            if RNA_result[0]["gene_name_ensembl"]!=".":
                                URL =URL+"," + "<a href=\"https://www.ncbi.nlm.nih.gov/gene/?term="+RNA_result[0]["gene_name_ensembl"]+"\">NCBI</a>"
                                URL =URL+"," + "<a href=\"https://ghr.nlm.nih.gov/search?query="+RNA_result[0]["gene_name_ensembl"]+"\">NIH - GHI</a>"
                            Final_result[Result_line].update({"External":URL})
                        else:
                            Final_result[Result_line].update({"gene_name_ensembl":".","gene_description":".","ensembl_gene_id":".","External":"."})
                            value=[]
                            value.extend(repeat(".",len(tissue_column_out.split(","))-1))
                            RNA_result=dict(zip(list(tissue_column_out.split(","))[1:],value))
                            Final_result[Result_line].update(RNA_result)



                        sql_command="SELECT %s FROM allele_frequency.%s WHERE (chr =\"%s\") AND (pos=\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % ("REVEL","REVEL",CHROM,record.POS,record.REF,record.ALT[idx])
                        con.execute(sql_command)
                        result = con.fetchall()

                        if len(result) > 0:
                            Final_result[Result_line].update(result[0])
                        else:
                            Final_result[Result_line].update({"REVEL":"."})




                        if len(request.form.getlist("ClinVar"))!=0:
                            sql_command = "SELECT `AlleleID`,`ClinicalSignificance`,`PhenotypeList` FROM allele_frequency.ClinVar WHERE (Chromosome =\"%s\") AND (Start=\"%s\") AND (ReferenceAllele = \"%s\") AND (AlternateAllele = \"%s\")" % (CHROM,record.POS,record.REF,record.ALT[idx])
                            con.execute(sql_command)
                            result = con.fetchall()
                            if len(result) > 0:
                                Final_result[Result_line].update(result[0])
                            else:
                                value = []
                                value.extend(repeat(".",3))
                                ClinVar_result = dict(zip(["AlleleID","ClinicalSignificance","PhenotypeList"],value))
                                Final_result[Result_line].update(ClinVar_result)
                        if len(request.form.getlist("CADD"))!=0:
                            Final_result[Result_line].update({"CADD_RawScore":".","CADD_PHRED":"."})


                        Result_line +=1
                if len(request.form.getlist("CADD"))!=0:
                    CADD_file.close()

            except Exception as e:
                #print(e)
                flash('File format not supported!')
                removeuserfile(User_id)
                if len(request.form.getlist("CADD"))!=0:
                    CADD_file.close()
                    return redirect(url_for('Variants_search'))
        else:
            flash('File format not supported!')
            removeuserfile(User_id)
            if len(request.form.getlist("CADD"))!=0:
                CADD_file.close()
            return redirect(url_for('Variants_search'))


       #ANNOVAR

        import subprocess
        import sys

        anno_data ="".join(list(map(lambda orig_string:","+orig_string,request.form.getlist("anno"))))
        anno_oper = ",f" * len(request.form.getlist("anno"))


        if file.filename == '':
            cmd='%s %s %s -buildver hg19 -remove --thread 2 -protocol ensGene,gerp++gt2%s -operation g,f%s -nastring . --outfile %s' % (ANNOVAR_path, os.path.join(app.config['UPLOAD_FOLDER'],User_id +".avinput"),humandb_path,anno_data,anno_oper, os.path.join(app.config['UPLOAD_FOLDER'],User_id))


        else:
            cmd='%s %s %s -buildver hg19 -remove --thread 2 -protocol ensGene,gerp++gt2%s -operation g,f%s -nastring . -vcfinput --outfile %s' % (ANNOVAR_path,os.path.join(app.config['UPLOAD_FOLDER'],filename),humandb_path,anno_data,anno_oper,os.path.join(app.config['UPLOAD_FOLDER'],User_id))

        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            #sys.exit(retcode)
            flash('File format not supported!')
            return redirect(url_for('Variants_search'))
        Current_line = 0

        with open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".hg19_multianno.txt"),"r") as annovar:
            for line in csv.DictReader(annovar,delimiter='\t'):
                annovar_result = [line[x] for x in ANNOVAR_column]
                Final_result[Current_line].update(dict(zip(ANNOVAR_column,annovar_result)))
                Current_line +=1


        if len(request.form.getlist("CADD"))!=0:

            cmd = [CADD_path,os.path.join(app.config['UPLOAD_FOLDER'],User_id+".vcf.gz"),os.path.join(app.config['UPLOAD_FOLDER'],User_id+".CADD.gz")]

            retcode = subprocess.call(cmd,env={"PATH":"/home/bioinfo/perl5/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/home/bioinfo/src/tabix/"})

            if retcode !=0: sys.exit(retcode)
            with gzip.open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".CADD.gz"), 'rb') as CADD_out:
                for line in CADD_out:
                    line = line.decode("utf-8").strip()
                    if line[0]!="#":
                        lines = line.split("\t")
                        line_dict_key = (lines[0],lines[1],lines[2],lines[3])
                        if line_dict_key in line_dict:
                            for result_line in line_dict[line_dict_key]:
                                Final_result[result_line].update({"CADD_RawScore":lines[4],"CADD_PHRED":lines[5]})


        if len(request.form.getlist("InterVar"))!=0:
            if file.filename == '':
                cmd = ["python2.7",InterVar_Path,"-b","hg19","-i",os.path.join(app.config['UPLOAD_FOLDER'],User_id +".avinput"),"-o",os.path.join(app.config['UPLOAD_FOLDER'],User_id)]
            else:
                cmd =["python2.7",InterVar_Path,"-b","hg19","-i",os.path.join(app.config['UPLOAD_FOLDER'],filename),"--input_type=VCF","-o",os.path.join(app.config['UPLOAD_FOLDER'],User_id)]

            retcode = subprocess.call(cmd,shell=True)

            if retcode != 0: sys.exit(retcode)
            Current_line = 0
            with open(os.path.join(app.config['UPLOAD_FOLDER'],User_id+".hg19_multianno.txt.intervar"),"r",encoding='utf-8',errors='ignore') as Intervar:
                for line in Intervar:
                    Intervar_type = ""
                    PS = 0
                    PM = 0
                    PP = 0
                    BS = 0
                    BP = 0
                    Current_line+=1
                    line = line.strip()
                    if Current_line!=1:
                        lines = line.split("\t")
                        intervar = lines[13].replace("InterVar: ","").replace(", ",",").split(" ")
                        if intervar[1] == "Likely" or intervar[1] == "Uncertain":
                            Intervar_type = " ".join(intervar[1:3])
                            PS = sum(list(map(int,intervar[4].split("=")[1][1:-1].split(","))))
                            PM = sum(list(map(int,intervar[5].split("=")[1][1:-1].split(","))))
                            PP = sum(list(map(int,intervar[6].split("=")[1][1:-1].split(","))))
                            BS = sum(list(map(int,intervar[8].split("=")[1][1:-1].split(","))))
                            BP = sum(list(map(int,intervar[9].split("=")[1][1:-1].split(","))))
                            if intervar[3].split("=")[1] == "1":
                                Intervar_type = Intervar_type + "(PVS1"
                            Intervar_type = Intervar_type + ",PS" + str(PS)
                            Intervar_type = Intervar_type + ",PM" + str(PM)
                            Intervar_type = Intervar_type + ",PP" + str(PP)
                            if intervar[7].split("=")[1] == "1":
                                Intervar_type = Intervar_type + "BA1"
                            Intervar_type = Intervar_type + ",BS" + str(BS)
                            Intervar_type = Intervar_type + ",BP" +str(BP)+")"



                        else:
                            Intervar_type = intervar[1]+"("
                            PS = sum(list(map(int,intervar[3].split("=")[1][1:-1].split(","))))
                            PM = sum(list(map(int,intervar[4].split("=")[1][1:-1].split(","))))
                            PP = sum(list(map(int,intervar[5].split("=")[1][1:-1].split(","))))
                            BS = sum(list(map(int,intervar[7].split("=")[1][1:-1].split(","))))
                            BP = sum(list(map(int,intervar[8].split("=")[1][1:-1].split(","))))
                            if intervar[2].split("=")[1] == "1":
                                Intervar_type = Intervar_type + "PVS1"
                            Intervar_type = Intervar_type + ",PS" + str(PS)
                            Intervar_type = Intervar_type + ",PM" + str(PM)
                            Intervar_type = Intervar_type + ",PP" + str(PP)
                            if intervar[6].split("=")[1] == "1":
                                Intervar_type = Intervar_type + "BA1"
                            Intervar_type = Intervar_type + ",BS" + str(BS)
                            Intervar_type = Intervar_type + ",BP" +str(BP)+")"

                        Final_result[Current_line-2].update({"InterVar":Intervar_type})


            db.close()
        """user table schema"""
        db=MySQLdb.connect(host="localhost",user="User_table",passwd="H58RrCeh69vtVuLq",db="user_table")
        con=db.cursor()

        if len(Final_result)!=0:
            Total_key = list(Final_result[0].keys())
            if "Index" not in Total_key:
                Total_key.append("Index")
            User_table_command = "CREATE TABLE `%s` (`id` int(11) NOT NULL auto_increment,%s,PRIMARY KEY  (`id`))" % (User_id,",".join(list(map(lambda orig:"`"+orig+"` TEXT NOT NULL",Total_key))))
            con.execute(User_table_command)
            User_insert = "INSERT INTO `%s` %s VALUES %s"

            """Add index"""
            for resaddID in Final_result:
                if "intergenic" in resaddID.values() or "intronic" in resaddID.values() or "synonymous SNV" in resaddID.values():
                    resaddID['Index'] = 0
                if resaddID['REVEL']!= ".":
                    if "nonsynonymous SNV" in resaddID.values() or float(resaddID['REVEL'])<REVEL_threshold:
                        resaddID['Index'] = 1
                    if "nonsynonymous SNV" in resaddID.values():
                        if float(resaddID['REVEL'])>REVEL_threshold or "splicing" in resaddID.values():
                            resaddID['Index'] = 2
                            if resaddID[target_tissue+"_FPKM"]!="." and resaddID[target_tissue+"_FPKM"]!=None:
                                if float(resaddID[target_tissue+"_FPKM"])>=0.5:
                                    resaddID['Index'] = 3
                            if resaddID['gerp++gt2']!=".":
                                resaddID['Index'] = 3
                elif resaddID['REVEL']== ".":
                    if "nonsynonymous SNV" in resaddID.values():
                        resaddID['Index'] = 1
                    if "splicing" in resaddID.values():
                        resaddID['Index'] = 2
                        if resaddID[target_tissue+"_FPKM"]!="." and resaddID[target_tissue+"_FPKM"]!=None:
                            if float(resaddID[target_tissue+"_FPKM"])>=0.5:
                                resaddID['Index'] = 3
                        if resaddID['gerp++gt2']!=".":
                                resaddID['Index'] = 3
                if "Index" not in resaddID:
                    resaddID['Index'] = 0
                column_key = str(tuple(resaddID.keys())).replace("\'","`")
                column_value = tuple(map(str,list(resaddID.values())))
                User_insert = "INSERT INTO `%s` %s VALUES %s" % (User_id,column_key,str(column_value))

                try:
                    con.execute(User_insert)
                    db.commit()
                #except:
                except Exception as e:
                    print(e)
                    db.rollback()
        else:
            db=MySQLdb.connect(host="localhost",user="User_table",passwd="H58RrCeh69vtVuLq",db="user_table")
            con=db.cursor()
            User_table_command = "CREATE TABLE `%s` (`id` int(11) NOT NULL auto_increment,%s,PRIMARY KEY  (`id`))" % (User_id,",".join(list(map(lambda orig:"`"+orig+"` TEXT NOT NULL",Final_result_title))))
            con.execute(User_table_command)
        User_table_command = "SELECT %s FROM `%s`" % (",".join(list(map(lambda orig:"`"+orig+"`",Final_result_title))),User_id)
        con.execute(User_table_command)
        Export_result = list(con.fetchall())
        Export_result.insert(0,Final_result_title)
        Current_time = time.strftime("%Y-%m-%dT%H:%M:%S")
        con.execute("INSERT INTO Table_history (User_id,Create_time,Now) VALUE (\"%s\",\"%s\",\"%s\")" % (User_id,Current_time,Current_time))



        db.close()





        """remove user files"""
        import glob
        for rmfile in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], User_id+"*")):
            os.remove(rmfile)




        return render_template("Teresult.html",Export_result=Export_result,results=Final_result,Gene_annotation=Gene_annotation,Population_allele_freq=Population_allele_freq,Predict=Predict,User_id=User_id,Clinical_interpretation=Clinical_interpretation)
    elif request.method == 'GET' and ("User_id" in request.values):
        print(request.values)
        User_id = request.values["User_id"]
        columns = request.values["table_column"].split(",")
        ajax_result = DataTablesServer(request,columns,User_id).outputResult()

        return ajax_result
    return render_template('Variants_search.html')

def removeuserfile(User_id):
    import glob
    for rmfile in glob.glob(os.path.join(app.config['UPLOAD_FOLDER'], User_id+"*")):
        os.remove(rmfile)
    return None

@app.route("/Tutorial")
def Tutorial():
    return render_template('Tutorial.html')

@app.route("/Contact")
def Contact():
    return render_template('Contact.html')
class DataTablesServer(object):
    def __init__( self, request, columns, collection):
        self.columns = columns
        self.index = "id"
        self.collection = collection
        self.request_values = request.values
        self.dbh = MySQLdb.connect(host="localhost",user="root",passwd="bioinfo",db="user_table")
        self.resultData = None
        self.cadinalityFiltered = 0
        self.cadinality = 0
        self.runQueries()
        self.outputResult()
    def outputResult(self):
        output = '{'
        output += '"draw": '+str(int(self.request_values['draw']))+', '
        output += '"recordsTotal": '+str(self.cardinality)+', '
        output += '"recordsFiltered": '+str(self.cadinalityFiltered)+', '
        output += '"data": [ '
        for row in self.resultData:
            output += '['
            for i in range( len(self.columns) ):
                output += '"'+row[ self.columns[i] ].replace('"','\\"')+'",'
            output = output[:-1]
            output += '],'
        output = output[:-1]
        output += '] }'
        return output

    def runQueries(self):
        dataCursor = self.dbh.cursor(cursorclass=MySQLdb.cursors.DictCursor)
        dataCursor.execute("""SELECT SQL_CALC_FOUND_ROWS %(columns)s FROM \
        `%(table)s` %(where)s %(order)s %(limit)s""" % \
        dict(columns=', '.join(list(map(lambda orig:"`"+orig+"`",self.columns))), table=self.collection ,where=self.filtering(),order=self.ordering(),limit=self.paging()))

        self.resultData = dataCursor.fetchall()
        cadinalityFilteredCursor = self.dbh.cursor()
        cadinalityFilteredCursor.execute( """ SELECT FOUND_ROWS() """ )
        self.cadinalityFiltered = cadinalityFilteredCursor.fetchone()[0]
        cadinalityCursor = self.dbh.cursor()
        cadinalityCursor.execute( """ SELECT COUNT(%s) FROM `%s` """ % (self.index, self.collection))
        self.cardinality = cadinalityCursor.fetchone()[0]
    def filtering(self):
        filter = ""
        if ( self.request_values.has_key('search[value]') ) and (self.request_values['search[value]'] != "" ):
            filter = "WHERE "
            for i in range(len(self.columns)):
                filter += "`%s` LIKE '%%%s%%' OR " % (self.columns[i],self.request_values['search[value]'])
            filter = filter[:-3]
        return filter
    def ordering( self ):
        order = ""
        if ( self.request_values['order[0][column]'] != "" ) and (int(self.request_values['order[0][column]']) > 0 ):

            order += "ORDER BY `%s` %s, " % (self.columns[int(self.request_values['order[0][column]'])],self.request_values['order[0][dir]'])
        return order[:-2]
    def paging(self):
        limit = ""
        if (self.request_values['start'] != "" ) and (self.request_values['length'] != -1 ):
            limit = "LIMIT %s, %s" % \
            (int(self.request_values['start']),int(self.request_values['length']))
        return limit



if __name__=='__main__':
    app.run(host='0.0.0.0',threaded=True)
