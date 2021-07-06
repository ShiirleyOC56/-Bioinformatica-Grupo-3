import os
from procesamientoImagen import operador_raiz,threholding, operador_exponencial, operador_logaritmico, histogram_equalization, constrast_streching, constrast_streching_out, power_raise
from alineamientoSecuencias import fastatoString,needleman_wunsch
from flask import Flask, render_template
from flask import url_for
from flask import redirect  
from flask import request
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = "./static/Fasta"

@app.route('/')
def home():
    return render_template('layout.html')  
    
@app.route('/formulario')
def upload_file():
    return render_template('layout.html')

@app.route("/raiz", methods=['POST'])
def endPointRaiz():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = operador_raiz(app.config['UPLOAD_FOLDER'], filename)

        return  render_template('layout.html', imagenR=imagen_resultado)

@app.route("/logaritmico", methods=['POST'])
def endPointLogaritmico():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = operador_logaritmico(app.config['UPLOAD_FOLDER'], filename)
 
        return  render_template('layout.html', imagenL=imagen_resultado)

@app.route("/exponencial", methods=['POST'])
def endPointExponencial():
    if request.method == 'POST':
        f = request.files['archivo']
        C = int(request.form['C'])

        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = operador_exponencial(app.config['UPLOAD_FOLDER'], filename, C)

        return  render_template('layout.html', imagenE=imagen_resultado)

@app.route("/histEquali", methods=['POST'])
def endPointHistEqua():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = histogram_equalization(app.config['UPLOAD_FOLDER'], filename)

        return  render_template('layout.html', imagenq=imagen_resultado)


@app.route("/contrastre", methods=['POST'])
def endPointContrastre():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = constrast_streching(app.config['UPLOAD_FOLDER'], filename)

        return  render_template('layout.html', imagenContraste=imagen_resultado)


@app.route("/contrasTreO", methods=['POST'])
def endPointContrasTreOu():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = constrast_streching_out(app.config['UPLOAD_FOLDER'], filename)

        return  render_template('layout.html', imagenendPct=imagen_resultado)


@app.route("/power", methods=['POST'])
def endPointRaisePower():
    if request.method == 'POST':
        f = request.files['archivo']
        c = int(request.form['b'])
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = power_raise(app.config['UPLOAD_FOLDER'], filename, c)

        return  render_template('layout.html', imagenRaisePower=imagen_resultado)


"""@app.route("/thresholding", methods=['POST'])
def endPointThresholding():
    if request.method == 'POST':
        f = request.files['archivo']
        filename = secure_filename(f.filename)
        full_filename = os.path.join(app.config['UPLOAD_FOLDER'], filename)
        f.save(full_filename)
        imagen_resultado = threholding(app.config['UPLOAD_FOLDER'], filename)

        return  render_template('layout.html', imagenThre=imagen_resultado)
"""
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

        score,alineamiento = needleman_wunsch(seq1,seq2,int(float(match)),int(float(mismatch)),int(float(gap)))
        return  render_template('layout.html', alineamientoGlobal=alineamiento,scoreM=score)

if __name__ == '__main__':
    app.run(debug=True)
