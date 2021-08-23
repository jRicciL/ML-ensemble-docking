import os
import pickle
import joblib
from functools import wraps

def run_or_load(func):
    '''Decorates a function with a "filename" parameter which is used to save a pickle
    file after run the function or directly load the file if it already exists'''
    def wrapper(filename, *args, **kwargs):
        if os.path.isfile(filename):
            print('File loaded:', filename)
            with open(filename, 'rb') as f:
                obj = pickle.load(f)
            return obj 
        else:
            obj = func(filename, *args, **kwargs)
            with open(filename, 'wb') as f:
                pickle.dump(obj, f)
            print('File saved:', filename)
            return obj
    return wrapper

def run_or_load_joblib(func):
    '''Decorates a function with a "filename" parameter which is used to save a pickle
    file after run the function or directly load the file if it already exists'''
    @wraps(func)
    def wrapper(filename, *args, **kwargs):
        if os.path.isfile(filename):
            print('File loaded:', filename)
            with open(filename, 'rb') as f:
                obj = joblib.load(f)
            return obj 
        else:
            obj = func(filename, *args, **kwargs)
            with open(filename, 'wb') as f:
                joblib.dump(obj, f)
            print('File saved:', filename)
            return obj
    return wrapper