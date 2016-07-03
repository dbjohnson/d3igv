import os
from flask import Flask
from flask import request
from flask import jsonify
from flask import send_from_directory


from app import read_generator


api = Flask(__name__)


DEFAULT_ARGS = {'read_error_rate': 0.03,
                'sequence_length': 80,
                'num_reads': 100,
                'max_SNPs': 4,
                'tumor_content': 0.8}


@api.route('/reads', methods=['POST'])
def get_reads():
    args = request.json.copy()
    for k, v in DEFAULT_ARGS.items():
        if k not in args:
            args[k] = v

    return jsonify(read_generator.generate_fake_reads(**args))


@api.route('/', methods=['GET'])
def index():
    return send_from_directory(os.path.dirname(__file__), 'index.html')
