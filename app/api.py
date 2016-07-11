import os
from flask import Flask
from flask import request
from flask import jsonify
from flask import send_from_directory


from app.bam import SimBAM


api = Flask(__name__)


DEFAULT_ARGS = {'read_error_rate': 0.003,
                'tumor_content': 0.3,
                'chrom': 'ch1',
                'start': 0,
                'end': 80,
                'num_reads': 100,
                'max_SNPs': 4}


@api.route('/reads', methods=['POST'])
def get_reads():
    args = request.json.copy()
    for k, v in DEFAULT_ARGS.items():
        if k not in args:
            args[k] = v

    simbam = SimBAM(args['read_error_rate'], args['tumor_content'])
    args.pop('read_error_rate')
    args.pop('tumor_content')
    return jsonify(simbam.fetch_segment(**args))


@api.route('/', methods=['GET'])
def default():
    return send_from_directory(os.path.dirname(__file__), 'index.html')


@api.route('/<path:filename>', methods=['GET'])
def static_pages(filename):
    return send_from_directory(os.path.dirname(__file__), filename)
