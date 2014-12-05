from copy import deepcopy
class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        #self.__dict__ = self


    def __getattr__(self, name):
        return self[name]

    def __setattr__(self, name, value):
        self[name] = value

    def __deepcopy__(self, memo):
        return AttrDict(deepcopy(dict(self), memo))