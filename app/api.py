import os
from flask import Flask
from flask import request
from flask import jsonify
from flask import send_from_directory


from app.bam import BAM
from app.bam import SimBAM


api = Flask(__name__)

simbam = SimBAM(read_error_rate=0.003, tumor_content=0.3)
bam = BAM('bam/demo.bam')


@api.route('/reads', methods=['POST'])
def get_reads():
    args = request.json
    sim = args.pop('sim', False)
    if sim:
        return jsonify(simbam.fetch_segment(**args))
    else:
        return jsonify(bam.fetch_segment(**args))


@api.route('/amplicons', methods=['GET'])
def load_bam():
    with open('bam/demo.bed', 'r') as fh:
        amplicons = [line.split() for line in fh]
    return jsonify(amplicons)


@api.route('/', methods=['GET'])
def default():
    return send_from_directory(os.path.dirname(__file__), 'index.html')


@api.route('/<path:filename>', methods=['GET'])
def static_pages(filename):
    return send_from_directory(os.path.dirname(__file__), filename)
