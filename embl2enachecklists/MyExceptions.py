#!/usr/bin/env python
'''
Class for custom exceptions
'''

###############
# AUTHOR INFO #
###############

__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>'
__copyright__ = 'Copyright (C) 2016-2017 Michael Gruenstaeudl'
__info__ = 'nex2embl'
__version__ = '2016.02.18.1100'

###########
# CLASSES #
###########


class MyException(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(1)

    def getErrorName(self):
        return "MyException"

class FileNotExist(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(2)

    def getErrorName(self):
        return "FileNotExist"

class WrongInputFile(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(3)

    def getErrorName(self):
        return "WrongInputFile"

class WrongOutputFile(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(4)

    def getErrorName(self):
        return "WrongOutputFile"

class ParserError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(5)

    def getErrorName(self):
        return "ParserError"

class ErrorNotFound(Exception):
    def __init__(self, value):
        self.value = "This Error does not exist and should not be seeable"

    def __str__(self):
        return(self.value)

    def getErrorNumber(self):
        return str(404)

    def getErrorName(self):
        return "ErrorNotFound"
