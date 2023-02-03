import hashlib
import json
import pickle
from datetime import datetime
import pandas as pd

def generate_uqid(thing) -> int:

    def json_default(thing):
        if isinstance(thing, datetime):
            return thing.isoformat(timespec='microseconds')

        if isinstance(thing, str):
            return pickle.dumps(thing)
        
        if isinstance(thing, int):
            return pickle.dumps(thing)

        raise TypeError(f"object of type {type(thing).__name__} not serializable")

    def json_dumps(thing):
        return json.dumps(
            thing,
            default=json_default,
            ensure_ascii=False,
            sort_keys=True,
            indent=None,
            separators=(',', ':'),
        )

    def get_hash(thing):
        return hashlib.md5(json_dumps(thing).encode('utf-8')).digest()

    return get_hash(thing).hex()