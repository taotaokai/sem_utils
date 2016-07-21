#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Convert XLS table to list

"""
import sys
from openpyxl import load_workbook

#====== parameters
xlsx_filename = str(sys.argv[1])

#====== read XLSX file
wb = load_workbook(filename=xlsx_filename, read_only=True)

for sheetname in wb.get_sheet_names():

    # get one worksheet (each regional network: Province)
    ws = wb[sheetname]

    # loop each line
    for row in ws.rows:
        if isinstance(row[0].value,int):
            name = row[1].value.encode('utf-8')
            country = str(row[2].value)
            subregion = str(row[7].value)
            tectonic = str(row[12].value)
            lat = float(row[8].value)
            lon = float(row[9].value)
            ele = float(row[10].value)
            print '%s|%s|%s|%s|%7.3f|%8.3f|%5.0f' % (
                    name,country,subregion,tectonic,lat,lon,ele)