import pandas as pd



def main(CSDF, CSThreshold):
    StackedCSDF = CSDF.stack()
    Qualified_Pairs = StackedCSDF[StackedCSDF < CSThreshold].to_dict()
    





















