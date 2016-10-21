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




UPLOAD_FOLDER='/home/bioinfo/Heart_gene_database/HeartInter/Uploads'
ALLOWED_EXTENSIONS=set(['txt','vcf'])


app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

"""
class fakefloat(float):
    def __init__(self, value):
        self._value = value
    def __repr__(self):
        return str(self._value)
for decimal data
def defaultencode(o):
    if isinstance(o, Decimal):
    # Subclass float with custom repr?
        return fakefloat(o)
    raise TypeError(repr(o) + " is not JSON serializable")"""
"""file reader"""
def allowed_file(filename):
    return '.' in filename and \
            filename.rsplit('.', 1)[1] in ALLOWED_EXTENSIONS

"""connect to the mySQL database
db=MySQLdb.connect(host="localhost",user="root",passwd="bioinfo",db="Heart_gene_expression")
con=db.cursor()
hsa_sql_command='''select
sample_info.description,sample_info.ref_ver,hsa_gene_info.gene_symbol,hsa_gene_info.Ensembl_gene_id,hsa_gene_info.description,hsa_gene_corr.chr,hsa_gene_corr.start,hsa_gene_corr.end,hsa_gene_corr.strand,hsa_gene_exp.FPKM,hsa_gene_exp.Level from sample_info,hsa_gene_info,hsa_gene_corr,hsa_gene_exp where sample_info.sample_id=hsa_gene_corr.ref_ver AND hsa_gene_corr.gene_id=hsa_gene_info.gene_id AND hsa_gene_exp.system_id=hsa_gene_info.gene_id '''

mus_sql_command='''select
    sample_info.description,sample_info.ref_ver,mus_gene_info.gene_symbol,mus_gene_info.Ensembl_gene_id,mus_gene_info.description,mus_gene_corr.Chr,mus_gene_corr.start,mus_gene_corr.end,mus_gene_corr.strand,mus_gene_exp.FPKM,mus_gene_exp.Level from sample_info,mus_gene_info,mus_gene_corr,mus_gene_exp where sample_info.sample_id=mus_gene_corr.ref_ver AND mus_gene_corr.gene_id=mus_gene_info.gene_id AND mus_gene_exp.system_id=mus_gene_info.gene_id '''

if (request.form.get('Zibrafish')):
    print('Zibrafish')

con.execute(hsa_sql_command)
hsa_records_dic=con.fetchall()

con.execute(mus_sql_command)
mus_records_dic=con.fetchall()

db.close()
"""
delimiters = '\t',';',"\r\n",',',' '
regexPattern = '|'.join(map(re.escape,delimiters))

def multiple_replace(string, rep_dict):
    pattern = re.compile("|".join([re.escape(k) for k in rep_dict.keys()]),re.M)
    return pattern.sub(lambda x: rep_dict[x.group(0)], string)

finalresults={}

@app.route("/")
def index():
    return render_template('index.html')

@app.route("/search")
def searchpage():
    return render_template('search.html')

@app.route("/search",methods=['GET','POST'])
def search():
    if request.method == "POST":
        db=MySQLdb.connect(host="localhost",user="user",passwd="bioinfo",db="Heart_gene_expression",cursorclass=MySQLdb.cursors.DictCursor)
        con=db.cursor()
        Search_text =multiple_replace(request.form['search_list'],{"\"":"","\'":"","`":"","%":""})
        Search_list = re.split(regexPattern,request.form['search_list'])

        FinalResult = {}

        if request.form.get('Human'):
            Ref_hum = "Ensembl_" + request.form.get('Human_ref') + "_RNA"
            #Ensembl gene information table columns
            RNA_field = request.form.getlist('RNA_field')
            protein_field = request.form.getlist('protein_field')


            #input gene symbol,gene symbol(ensembl)>gene symbol(NCBI)>aliases
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
                    FinalResult[Ensembl_Res["gene_name(ensembl)"]] = Ensembl_Res
                    Ensembl_Set.add(Ensembl_Res["gene_name(ensembl)"])

                if len(set(Search_list)-Ensembl_Set)!=0:
                    in_p = ', '.join(list(map(lambda x: '\'' + x + '\'',set(Search_list)-Ensembl_Set)))
                    NCBI_Gene = 'SELECT * FROM `%s` WHERE `gene_name(NCBI)` IN (%s)'
                    NCBI_Gene = NCBI_Gene % (Ref_hum,in_p)
                    con.execute(NCBI_Gene)
                    Result = con.fetchall()
                    for NCBI_Res in Result:
                        FinalResult[NCBI_Res["gene_name(NCBI)"]] = NCBI_Res
                        NCBI_Set.add(NCBI_Res["gene_name(NCBI)"])
                aliase = set(Search_list) - Ensembl_Set - NCBI_Set
                if len(set(aliase))!=0:
                    for symbol in aliase:
                        aliase_Gene = 'SELECT * FROM `%s` WHERE FIND_IN_SET("%s",REPLACE(`aliases`,"|",","))' % (Ref_hum,symbol)
                        con.execute(aliase_Gene)
                        Result = con.fetchall()
                        FinalResult[symbol] = Result
                        aliase_Set.add(symbol)
                NoData = aliase - aliase_Set
                if len(NoData)!=0:
                    for no in NoData:
                        FinalResult[no]=None

                print(FinalResult)



    return render_template('search.html')


@app.route('/vcf', methods=['GET','POST'])
def upload_file():
    if request.method == 'POST':
        file = request.files['file']
        if file and allowed_file(file.filename):
            filename = secure_filename(file.filename)
            file.save(os.path.join(app.config['UPLOAD_FOLDER'], filename))
            vcf_reader = vcf.Reader(open(filename, 'r'))
            for record in vcf_reader:
                print(record)

            return redirect(url_for('uploaded_file',filename=filename))
    return render_template('vcf.html')

@app.route('/uploads/<filename>')
def uploaded_file(filename):
    return send_from_directory(app.config['UPLOAD_FOLDER'],filename)
"""
@app.route('/variant',methods=['GET','POST'])
def variant():
"""




if __name__=='__main__':
    app.run(host='0.0.0.0',port=5000,debug=True)

#remember "debug=True" should be deleted!!
