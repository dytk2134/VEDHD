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
    freq_result = []
    sql_command="SELECT %s FROM allele_frequency.%s WHERE (chr =\"%s\") AND (pos=\"%s\") AND (ref = \"%s\") AND (alt = \"%s\")" % (column,table,chr,pos,ref,alt)
    con.execute(sql_command)
    result = con.fetchall()
    if len(result) > 0:
        freq_result.extend(list(result[0]))
    else:
        from itertools import repeat
        freq_result.extend(repeat(".",poplen))
    return freq_result


@app.route("/")
def index():
    return render_template('index.html')

@app.route("/search",methods=['GET','POST'])
def search():
    if request.method == "POST":
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Heart_gene_expression")
        con=db.cursor()
        Search_text=multiple_replace(request.form['search_list'],{"\"":"","\'":"","`":"","%":""})
        Search_text=Search_text.upper()
        Search_list =list(filter(None,re.split(regexPattern,Search_text)))
        FinalResult = []
        if request.form.get('Human'):
            FinalResult_num = 0
            FinalResult_title = ["user_input","gene_name(ensembl)","gene_name(NCBI)","aliases","Ensembl_id_GRCh38.p5","Ensembl_id_GRCh37","chr","gene_start","gene_end","strand","gene_description"]
            if request.form["output"] == "RNA":
                FinalResult_title.extend(["tissue","FPKM","RNA_level"])
            else:
                FinalResult_title.extend(["tissue","cell_type","Expression_type","Protein_Level","Reliability"])
            Ref_hum = "Ensembl_" + request.form.get('Human_ref') + "_" + request.form["output"]
            if request.form.get('SearchType') == "Symbol":
                Ensembl_Set = set()
                NCBI_Set = set()
                aliase_Set = set()
                Ensembl_Gene = 'SELECT * FROM `%s` WHERE `gene_name(ensembl)` IN (%s)'
                in_p = ', '.join(list(map(lambda x: '\'' + x + '\'', Search_list)))
                Ensembl_Gene = Ensembl_Gene % (Ref_hum,in_p)
                con.execute(Ensembl_Gene)
                Result = con.fetchall()

                for Ensembl_Res in Result:
                    FinalResult.append([Ensembl_Res[0]])
                    FinalResult[FinalResult_num].extend(list(Ensembl_Res))
                    Ensembl_Set.add(Ensembl_Res[0])
                    FinalResult_num +=1

                if len(set(Search_list)-Ensembl_Set)!=0:
                    in_p = ', '.join(list(map(lambda x: '\'' + x +'\'',set(Search_list)-Ensembl_Set)))
                    NCBI_Gene = 'SELECT * FROM `%s` WHERE `gene_name(NCBI)` IN (%s)'
                    NCBI_Gene = NCBI_Gene % (Ref_hum,in_p)
                    con.execute(NCBI_Gene)
                    Result = con.fetchall()
                    for NCBI_Res in Result:
                        FinalResult.append([NCBI_Res[1]])
                        FinalResult[FinalResult_num].extend(list(NCBI_Res))
                        NCBI_Set.add(NCBI_Res[0])
                        FinalResult_num +=1

                aliase = set(Search_list) - Ensembl_Set - NCBI_Set
                if len(set(aliase))!=0:
                    for symbol in aliase:
                        aliase_Gene = 'SELECT * FROM `%s` WHERE FIND_IN_SET("%s",REPLACE(`aliases`,"|",","))' % (Ref_hum,symbol)
                        con.execute(aliase_Gene)
                        Result = con.fetchall()
                        if len(Result) != 0:
                            FinalResult.append([symbol])
                            if len(Result)>1:
                                for aliase in Result:
                                    FinalResult[FinalResult_num].extend(list(aliase))
                                    FinalResult_num +=1
                                aliase_Set.add(sumbol)
                            elif len(Result) == 1:
                                FinalResult[FinalResult_num].extend(list(Result[0]))
                                aliase_Set.add(symbol)
                                FinalResult_num +=1
                NoData = aliase - aliase_Set
                if len(NoData)!=0:
                    from itertools import repeat
                    for no in NoData:
                        FinalResult.append([no])
                        FinalResult[FinalResult_num].extend(repeat("-",len(FinalResult_title)-1))
                        FinalResult_num +=1

            elif request.form.get('SearchType') == "Ensembl_id":
                Ensembl_Set = set()
                if request.form.get('Human_ref') == "GRCh38":
                    Ensembl_id_ver = "GRCh38.p5"
                else:
                    Ensembl_id_ver = "GRCh37"
                sql_command = 'SELECT * FROM `%s` WHERE `%s` in (%s)'
                in_p = ', '.join(list(map(lambda x: '\'' + x + '\'',Search_list)))
                sql_command = sql_command % (Ref_hum,Ensembl_id_ver,in_p)
                con.execute(sql_command)
                Result = con.fetchall()
                for Ensembl_Res in Result:
                    if Ensembl_id_ver == "GRCh38.p5":
                        FinalResult.append([Ensembl_Res[3]])
                        Ensembl_Set.add(Ensembl_Res[3])
                    else:
                        FinalResult.append([Ensembl_Res[4]])
                        Ensembl_Set.add(Ensembl_Res[4])
                    FinalResult[FinalResult_num].extend(list(Ensembl_Res))
                    FinalResult_num +=1
                if len(set(Search_list)-Ensembl_Set)!=0:
                    NoData = set(Search_list) - Ensembl_Set
                    from itertools import repeat
                    for no in NoData:
                        FinalResult.append([no])
                        FinalResult[FinalResult_num].extend(repeat("-",len(FinalResult_title)-1))
                        FinalResult_num +=1
            db.close()
        return render_template("results.html",results=FinalResult,keys=FinalResult_title)
    return render_template('new_search.html')


@app.route('/vcf', methods=['GET','POST'])
def upload_file():
    if request.method == 'POST':

        """Connect to MySQL database"""
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo")
        con=db.cursor()
        """user input information"""
        Output_format = request.form["format"]
        Genomes_population = request.form.getlist("Genomes")#1000 Genomes population list
        JPN_population = request.form.getlist("JPN")
        ESP_population = request.form.getlist("ESP")

        Final_result = [] #存最後的結果
        Final_result_title = ["chr","pos","ref","alt"]
        if len(Genomes_population)!=0:
            Genomes_population_column = ",".join(list(map(lambda orig_string:orig_string + "_Ref_" + Output_format +","+orig_string + "_Alt_"+ Output_format,Genomes_population)))
            Final_result_title.extend(Genomes_population_column.split(","))
        if "1KJPN" in JPN_population:
            JPN_1_column = "1KJPN_Ref_" + Output_format + "," + "1KJPN_Alt_" +Output_format
            Final_result_title.extend(JPN_1_column.split(","))
        if "2KJPN" in JPN_population:
            JPN_2_column = "2KJPN_Ref_" + Output_format + "," + "2KJPN_Alt_" +Output_format
            Final_result_title.extend(JPN_2_column.split(","))
        if len(ESP_population) !=0:
            ESP_population_column = ",".join(list(map(lambda orig_string:orig_string + "_Ref_" +Output_format +","+orig_string + "_Alt_" +Output_format,ESP_population)))
            Final_result_title.extend(ESP_population_column.split(","))
        if len(request.form.getlist("REVEL")) != 0:
            Final_result_title.extend(["aaref","aaalt","REVEL"])
        Final_result_title.extend(["gene_name","description","FPKM(Heart muscle)","Ranking/Total"])

        Result_line = 0


        #if user input file
        file = request.files['file']
        #如果使用者上傳檔案，則使用此檔案的資訊
        if file.filename == '':
            variant_text = multiple_replace(request.form['variants_list'],{"\"":"","\'":"","`":"","%":""})
            variant_text = variant_text.upper()
            variant_list =list(filter(None,re.split(regexPattern,variant_text)))
            vartext_file = open(os.path.join(app.config['UPLOAD_FOLDER'],"user_input.avinput"),"w")

            for var in variant_list:
                vars = re.compile("([\dXYM]+):(\d+)([ATCG]+)>([ATCG]+)").split(var)
                if len(vars)<6:
                    Final_result.append([var[0],".",".","."])
                    vartext_file.write(var[0]+"\t.\t.\t.\n")
                else:
                    Final_result.append(vars[1:5])
                    vartext_file.write("chr"+vars[1]+"\t"+vars[2]+"\t"+str(int(vars[2])+len(vars[3])-1)+"\t"+vars[3]+"\t"+vars[4]+"\n")


                if len(Genomes_population)!=0:
                    Final_result[Result_line].extend(get_freq(con,"1000Genomes_5pop_"+Output_format,Genomes_population_column,len(Genomes_population)*2,vars[1],vars[2],vars[3],vars[4]))
                if "1KJPN" in JPN_population:
                    Final_result[Result_line].extend(get_freq(con,"1KJPN_"+Output_format,JPN_1_column,2,vars[1],vars[2],vars[3],vars[4]))
                if "2KJPN" in JPN_population:
                    Final_result[Result_line].extend(get_freq(con,"2KJPN_"+Output_format,JPN_2_column,2,vars[1],vars[2],vars[3],vars[4]))
                if len(ESP_population) !=0:
                    Final_result[Result_line].extend(get_freq(con,"ESP_"+Output_format,ESP_population_column,len(ESP_population)*2,vars[1],vars[2],vars[3],vars[4]))

                if len(request.form.getlist("REVEL"))!=0:
                    Final_result[Result_line].extend(get_freq(con,"REVEL","aaref,aaalt,REVEL",3,vars[1],vars[2],vars[3],vars[4]))
                sql_command = "SELECT `gene_name`, `description`,`FPKM(heart muscle)`,`Ranking/Total` FROM Heart_gene_expression.Proteinaltas_RNA WHERE chr = \"%s\" AND %s BETWEEN Proteinaltas_RNA.start AND Proteinaltas_RNA.end" % (vars[1],vars[2])
                print(sql_command)
                con.execute(sql_command)
                RNA_result = con.fetchall()
                if len(RNA_result) > 0:
                    Final_result[Result_line].extend(list(RNA_result[0]))
                else:
                    from itertools import repeat
                    Final_result[Result_line].extend(repeat(".",4))
                Result_line +=1
            vartext_file.close()

        elif file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'],filename))
            vcf_reader = vcf.Reader(open(os.path.join(app.config['UPLOAD_FOLDER'],filename)))
            for record in vcf_reader:
                for idx in range(0,len(record.ALT)):
                    Final_result.append([record.CHROM,record.POS,record.REF,record.ALT[idx]])
                    if len(Genomes_population)!=0:
                        Final_result[Result_line].extend(get_freq(con,"1000Genomes_5pop_"+Output_format,Genomes_population_column,len(Genomes_population)*2,record.CHROM,record.POS,record.REF,record.ALT[idx]))

                    if "1KJPN" in JPN_population:
                        Final_result[Result_line].extend(get_freq(con,"1KJPN_"+Output_format,JPN_1_column,2,record.CHROM,record.POS,record.REF,record.ALT[idx]))
                    if "2KJPN" in JPN_population:
                       Final_result[Result_line].extend(get_freq(con,"2KJPN_"+Output_format,JPN_2_column,2,record.CHROM,record.POS,record.REF,record.ALT[idx]))
                    if len(ESP_population) !=0:
                        Final_result[Result_line].extend(get_freq(con,"ESP_"+Output_format,ESP_population_column,len(ESP_population)*2,record.CHROM,record.POS,record.REF,record.ALT[idx]))
                    if len(request.form.getlist("REVEL"))!=0:
                        Final_result[Result_line].extend(get_freq(con,"REVEL","aaref,aaalt,REVEL",3,record.CHROM,record.POS,record.REF,record.ALT[idx]))
                       #RNA expression
                    sql_command = "SELECT `gene_name`, `description`,`FPKM(heart muscle)`,`Ranking/Total` FROM Heart_gene_expression.Proteinaltas_RNA WHERE chr = \"%s\" AND %s BETWEEN Proteinaltas_RNA.start AND Proteinaltas_RNA.end" % (record.CHROM,record.POS)
                    con.execute(sql_command)
                    RNA_result = con.fetchall()
                    if len(RNA_result) > 0:
                        Final_result[Result_line].extend(list(RNA_result[0]))
                    else:
                        from itertools import repeat
                        Final_result[Result_line].extend(repeat(".",4))

                    Result_line +=1
        print(Final_result)

#ANNOVAR
        import subprocess
        if file.filename == '':
            cmd='/home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/table_annovar.pl %s /home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/humandb/ -buildver hg19 -remove -protocol refGene,gerp++gt2 -operation g,f -nastring . --outfile %s' % (os.path.join(app.config['UPLOAD_FOLDER'],"user_input.avinput"),os.path.join(app.config['UPLOAD_FOLDER'],"output"))
            print(cmd)
        else:
            cmd='/home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/table_annovar.pl %s /home/bioinfo/Heart_gene_database/HeartInter/code/tool/annovar/humandb/ -buildver hg19 -remove -protocol refGene,gerp++gt2 -operation g,f -nastring . -vcfinput --outfile %s' % (os.path.join(app.config['UPLOAD_FOLDER'],filename),os.path.join(app.config['UPLOAD_FOLDER'],"output"))

        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0: sys.exit(retcode)
        Current_line = 0
        Final_result_title.extend(["Func.refGene","ExonicFunc.refGene","AAChange.refGene","gerp++gt2"])
        with open(os.path.join(app.config['UPLOAD_FOLDER'],"output.hg19_multianno.txt"),"r") as annovar:
            for line in annovar:
                Current_line+=1
                line = line.strip()
                if Current_line !=1:
                    lines = line.split("\t")
                    Final_result[Current_line-2].extend([lines[5],lines[8],lines[9],lines[10]])

        db.close()
            #os.remove(os.path.join(app.config['UPLOADED_ITEMS_DEST'], filename))
        return render_template("Teresult.html",results = Final_result,keys = Final_result_title)
    return render_template('vcf.html')


if __name__=='__main__':
    app.run(host='0.0.0.0',port=5000,debug=True)
