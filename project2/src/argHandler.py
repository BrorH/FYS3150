
"""
Easily handles arguments passed via terminal. As long as arguemnts are numbers and variable have names this is perfect
balsadhaskdhkasd
fhsdfsrf
fill meeee

"""
class argHandler:
    def __init__(self, defaults):
        self.args = {}
        self.defaults = {}
        for key in defaults.keys():
            # makes all self.defaults lists. neccasairy for later work
            if not isinstance(defaults[key], list):
                self.defaults[key] = [defaults[key]]
    

    def splitrange(self, string):
        # turns a string-formated range, i.e 10:101:10 into a python list list(range(10,101, 10))
        # can be passed as 10:101 (interval as 1) or 10:101:23 (interval as 23)
        splitted = string.split(":")
        if len(splitted) == 2:
            return self.splitrange(string+":1")
        
        return list(range(int(splitted[0]), int(splitted[1]), int(splitted[2])))
    

    def splitlist(self, string):
        # splits a string of comma separated values into python list
        return [float(obj) for obj in string.split(",")]

    def parse(self, args):
        self.args = {}
        for arg in args[1:]:
            for argname in self.defaults.keys():
                eqIdx = arg.index("=")
                if arg[:eqIdx] == argname:
                    if "," in arg[eqIdx+1:]:
                        self.args[argname] = self.splitlist(arg[eqIdx+1:])
                    elif ":" in arg[eqIdx+1:]:
                        self.args[argname] = self.splitrange(arg[eqIdx+1:])
                    else:
                        self.args[argname] = [arg[eqIdx+1:]]
        for undefined in list(set(self.defaults.keys())- set(self.args.keys())):
            print(f"{undefined} defaulted to {self.defaults[undefined][0]}")
            self.args[undefined] = self.defaults[undefined]
        return self.args