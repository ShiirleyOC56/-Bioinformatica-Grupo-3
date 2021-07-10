from itertools import product
import os, glob
from alineamientoSecuencias import fastatoString,needleman_wunsch1,sw,bdtoarray,BLAST,MUSCLE,jukes_cantor,K2Pdistance,alpha_labels,UPGMA
from alineamientoSecuencias import TNdistance, Tamuradistance
from flask import Flask, render_template
from flask import url_for
from flask import redirect  
from flask import request
from werkzeug.utils import secure_filename
os.environ['QT_QPA_PLATFORM']='offscreen'
from ete3 import Tree
from IPython.display import display, Image

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
        secuencia2 = request.files['archivo2']
        stri1=request.form.get('entrada1')
        stri2=request.form.get('entrada2')
        match = request.form['match']
        mismatch = request.form['mismatch']
        gap = request.form['gapextend']
        if (secuencia1 and secuencia2):
            filename1 = secure_filename(secuencia1.filename)
            full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
            secuencia1.save(full_filename)

            filename2 = secure_filename(secuencia2.filename)
            full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
            secuencia2.save(full_filename2)

            seq1 = fastatoString("./static/Fasta/"+filename1)
            seq2 = fastatoString("./static/Fasta/"+filename2)
            
            print("SI EXISTEN LOS ARCHIVOS")
            alineamiento = needleman_wunsch1(seq1,seq2,int(float(match)),int(float(mismatch)),int(float(gap)))
            return  render_template('layout.html', tupla=alineamiento)
        elif (stri1 and stri2):
            print(stri1)
            print(stri2)
            alineamiento = needleman_wunsch1(stri1,stri2,int(float(match)),int(float(mismatch)),int(float(gap)))
            return  render_template('layout.html', tupla=alineamiento)



@app.route("/AlineamientoLocal", methods=['POST'])
def endPointAlineamientoLocal():
    if request.method == 'POST':        
        secuencia1 = request.files['archivo3']
        secuencia2 = request.files['archivo4']
        stri3 = request.form.get('entrada3')
        stri4= request.form.get('entrada4')
        gap = request.form['gapextend2']
        match = request.form['match2']
        mismatch = request.form['mismatch2']

        if (secuencia1 and secuencia2):
            filename1 = secure_filename(secuencia1.filename)
            full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
            secuencia1.save(full_filename)

            filename2 = secure_filename(secuencia2.filename)
            full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
            secuencia2.save(full_filename2)

            seq1 = fastatoString("./static/Fasta/"+filename1)
            seq2 = fastatoString("./static/Fasta/"+filename2)

            #alineamiento es un diccionario
            alineamientoL = sw(seq1,seq2,int(float(gap)),int(float(match)),int(float(mismatch)))
            return  render_template('layout.html', dic=alineamientoL)
        elif (stri3 and stri4):
            alineamientoL = sw(stri3,stri4,int(float(match)),int(float(mismatch)),int(float(gap)))
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

@app.route("/TamuraModel", methods=['POST'])
def endPointTamuraModel():
    if request.method == 'POST':
        secuencia1 = request.files['archivoAA']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivoBB']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        secuencia2.save(full_filename2)

        seq1 = fastatoString("./static/Fasta/"+filename1)
        seq2 = fastatoString("./static/Fasta/"+filename2)
        
        tamuraD = Tamuradistance(seq1, seq2)
        return  render_template('layout.html', distanceT=tamuraD)

@app.route("/TajimaModel", methods=['POST'])
def endPointTajimaModel():
    if request.method == 'POST':
        secuencia1 = request.files['archivoAAA']
        filename1 = secure_filename(secuencia1.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename1)
        secuencia1.save(full_filename)

        secuencia2 = request.files['archivoBBB']
        filename2 = secure_filename(secuencia2.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        secuencia2.save(full_filename2)

        seq1 = fastatoString("./static/Fasta/"+filename1)
        seq2 = fastatoString("./static/Fasta/"+filename2)
        
        tajimaD = TNdistance(seq1, seq2)
        print("tajimaaaaaaaaaaaaaaaa",tajimaD)
        return  render_template('layout.html', distanceTj=tajimaD)


@app.route("/UPGMA", methods=['POST'])
def endPointUPGMA():
    if request.method == 'POST':
        M = []
        sublista1 = []
        """a1= request.form['parametro1']
        a2= request.form.get('a2')
        a3= request.form.get('a3')
        a4= request.form.get('a4')
        sublista1.append(int(float(alfa1)))
        sublista1.append(int(float(a2)))
        sublista1.append(int(float(a3)))
        sublista1.append(int(float(a4)))"""
        M.append(sublista1)

        sublista2 = []
        b1= request.form['parametro5']
        """b2= request.form.get('b2')
        b3= request.form.get('b3')
        b4= request.form.get('b4')"""
        sublista2.append(int(float(b1)))
        """sublista2.append(int(float(b2)))
        sublista2.append(int(float(b3)))
        sublista2.append(int(float(b4)))"""
        M.append(sublista2)

        sublista3 = []
        c1= request.form['parametro9']
        c2= request.form['parametro10']
        """c3= request.form.get('c3')
        c4= request.form.get('c4')"""
        sublista3.append(int(float(c1)))
        sublista3.append(int(float(c2)))
        """sublista3.append(int(float(c3)))
        sublista3.append(int(float(c4)))"""
        M.append(sublista3)

        sublista4 = []
        d1= request.form['parametro13']
        d2= request.form['parametro14']
        d3= request.form['parametro15']
        """d4= request.form.get('d4')"""
        sublista4.append(int(float(d1)))
        sublista4.append(int(float(d2)))
        sublista4.append(int(float(d3)))
        """sublista4.append(int(float(d4)))"""
        M.append(sublista4)

        M_labels = alpha_labels("A", "D")

        upgma = UPGMA(M, M_labels)
        strr=upgma + ";"
        t = Tree(strr)# aqui copia lo que te salio en el upgma
        #t.render("./static/ArbolesFilogeneticos/figurita.png")
        t.render("./static/ArbolesFilogeneticos/figurita.png")
        img = "./static/ArbolesFilogeneticos/figurita.png"
        print(img)
        return  render_template('layout.html', up=upgma,img=img)

@app.route("/NeighborJoining", methods=['POST'])
def endPointNeighborJoining():
    if request.method == 'POST':
        atxt = request.files['archivotxt']
        filename2 = secure_filename(atxt.filename)
        full_filename2 = os.path.join(app.config['UPLOAD_FOLDER'], filename2)
        atxt.save(full_filename2)

        nnParametro = request.form['nn']
        llParametro = request.form['ll']

        dis_data = open(filename2)

        dis_map = {}
        i = 0
        for line in dis_data:
            split_line = line.split()
            print (split_line)
            for j in range(len(split_line)):
                dis_map[(i,j)] = int(split_line[j])
            i+=1
        dis_data.close()

        nj = neighborjoining(dis_map, int(float(nn)),int(float(ll)))

        return  render_template('layout.html', numfc=filascolumnas)

if __name__ == '__main__':
    app.run(debug=True)
