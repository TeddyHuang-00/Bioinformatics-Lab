import streamlit as slt

# import requests
# from bs4 import BeautifulSoup

# baseUrl = "http://gepia.cancer-pku.cn/detail.php"
# query = {"gene": "ERBB2", "tag": "expdiy"}

# res = requests.post(url=baseUrl, data=query)
# # response.text
# # slt.markdown(response.text[17:], unsafe_allow_html=True)
# # response.text[17:]
# soup = BeautifulSoup(res.text, "html.parser")
# soup.current_data
# with open("temp.html", "wb") as fout:
#     fout.write(res.content)

import gepia

bp = gepia.boxplot()
slt.text(str(bp.showParams()))

slt.write(gepia.CANCERType)
bp.setParam("dataset", gepia.CANCERType)
result = bp.query()

# IFrame(result, width=500, height=500)
