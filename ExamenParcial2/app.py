from itertools import product
import os, glob
from procesamientoImagen import operador_raiz,threholding, operador_exponencial, operador_logaritmico, histogram_equalization, constrast_streching, constrast_streching_out, power_raise
from alineamientoSecuencias import fastatoString,needleman_wunsch,sw,bdtoarray,BLAST,MUSCLE,jukes_cantor,K2Pdistance
from flask import Flask, render_template
from flask import url_for
from flask import redirect  
from flask import request
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = "./static/Fasta"
app.config['UPLOAD_FOLDER1'] = "./static/bdFasta"
@app.route('/')
def home():
    return render_template('layout.html')  
    
@app.route('/formulario')
def upload_file():
    return render_template('layout.html')


@app.route("/AlineamientoGlobal", methods=['POST'])
def endPointAlineamientoGlobal():
    if request.method == 'POST':
        secuencia1 = request.files['archivo1']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivo2']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        secuencia2.save(full_filename2)

        seq1 = fastatoString("./static/Fasta/"+filename1)
        seq2 = fastatoString("./static/Fasta/"+filename2)
        gap = request.form['gapextend']
        match = request.form['match']
        mismatch = request.form['mismatch']
        #alineamiento es un diccionario
        alineamiento = needleman_wunsch(seq1,seq2,int(float(match)),int(float(mismatch)),int(float(gap)))
        return  render_template('layout.html', dic=alineamiento)

@app.route("/AlineamientoLocal", methods=['POST'])
def endPointAlineamientoLocal():
    if request.method == 'POST':
        secuencia1 = request.files['archivo3']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivo4']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        secuencia2.save(full_filename2)

        seq1 = fastatoString("./static/Fasta/"+filename1)
        seq2 = fastatoString("./static/Fasta/"+filename2)
        gap = request.form['gapextend2']
        match = request.form['match2']
        mismatch = request.form['mismatch2']
        #alineamiento es un diccionario
        alineamientoL = sw(seq1,seq2,int(float(gap)),int(float(match)),int(float(mismatch)))
        return  render_template('layout.html', dic=alineamientoL)

@app.route("/BLAST", methods=['POST'])
def endPointBLAST():
    if request.method == 'POST':
        secuencia1 = request.files['archivo5']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER1'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivo6']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER1'], filename2)
        secuencia2.save(full_filename2)

        secuencia3 = request.files['archivo7']
        filename3 = secure_filename(secuencia3.filename)
        full_filename3 = os.path.join(app.config['UPLOAD_FOLDER1'], filename3)
        secuencia3.save(full_filename3)

        secuencia4 = request.files['archivo8']
        filename4 = secure_filename(secuencia4.filename)
        full_filename4 = os.path.join(app.config['UPLOAD_FOLDER1'], filename4)
        secuencia4.save(full_filename4)

        secuencia5 = request.files['archivo9']
        filename5 = secure_filename(secuencia5.filename)
        full_filename5 = os.path.join(app.config['UPLOAD_FOLDER1'], filename5)
        secuencia5.save(full_filename5)

        #SECUENCIA DE PRUEBA
        secuencia6 = request.files['archivo10']
        filename6 = secure_filename(secuencia6.filename)
        full_filename6 = os.path.join(app.config['UPLOAD_FOLDER'], filename6)
        secuencia6.save(full_filename6)

        bd = bdtoarray("./static/bdFasta/")
        query = fastatoString("./static/Fasta/"+filename6)

        blast = BLAST(query,bd)
        return  render_template('layout.html', dic=blast)

@app.route("/Muscle", methods=['POST'])
def endPointMuscle():
    if request.method == 'POST':
        secuencia1 = request.files['archivo11']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        alineamientoL = sw(seq1,seq2,int(float(gap)),int(float(match)),int(float(mismatch)))
        return  render_template('layout.html', dic=alineamientoL)

@app.route("/JukesCantor", methods=['POST'])
def endPointJukesCantor():
    if request.method == 'POST':
        secuencia1 = request.files['archivo12']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)
        dir = "./static/Fasta/"+filename1
        jukes = jukes_cantor(dir)
        return  render_template('layout.html', arr=jukes)

@app.route("/KimuraModel", methods=['POST'])
def endPointKimuraModel():
    if request.method == 'POST':
        secuencia1 = request.files['archivoA']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivoB']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        secuencia2.save(full_filename2)

        seq1 = fastatoString("./static/Fasta/"+filename1)
        seq2 = fastatoString("./static/Fasta/"+filename2)
        
        kimuraD = K2Pdistance(seq1, seq2)
        return  render_template('layout.html', distance=kimuraD)

if __name__ == '__main__':
    app.run(debug=True)
