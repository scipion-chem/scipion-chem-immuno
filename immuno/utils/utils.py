# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import time, os, requests
from Bio import SeqIO

from ..constants import EVAL_PARAM_MAP

def runEpitopeSelection(softwareName, argsDic, browserData={}):
  ''' Run an epitope selector program with the specified arguments and parse the results
  :param softwareName: Selector software to call
  :param argsDic: dictionary containing the arguments for the selector. Keys must be the ones expected by the program
  :return: {seq_id: {(position, epitopeString): meanScore}}
  '''
  if softwareName.lower() == 'abcpred':
    protsDic = parseInputProteins(argsDic['i'])
    if browserData:
      epiDic = callABCpredSelenium(protsDic, browserData, argsDic)
    else:
      epiDic = callABCpred(protsDic, argsDic)

  elif softwareName.lower() == 'lbtope':
    protsDic = parseInputProteins(argsDic['i'])
    epiDic = callLBtope(protsDic, browserData, argsDic)

  return epiDic



def parseInputProteins(faFile):
  '''Uses BioPython to parse a fasta file and return it as dictionary
  :param faFile: input fasta filename
  :return: {seqName1: seqStr1, ...}
  '''
  faDic = {}
  with open(faFile) as handle:
    for values in SeqIO.FastaIO.SimpleFastaParser(handle):
        faDic[values[0]] = values[1]
  return faDic


def reportPoolStatus(poolDic):
  '''Check the status of the AsynPool objects stored as values of the dictionary and reports when they finish
  '''
  ready = []
  while len(ready) < len(poolDic):
    time.sleep(5)
    for evalSoft, po in poolDic.items():
      if po.ready() and evalSoft not in ready:
        ready.append(evalSoft)
        print(f'{evalSoft} execution finished ({len(ready)} / {len(poolDic)})')



def divide_chunks(iter, chunkSize):
  '''Divides an iterable into chunks of size chunkSize'''
  chunks = []
  for i in range(0, len(iter), chunkSize):
    chunks.append(iter[i:i + chunkSize])
  return chunks

def buildSeqFasta(seqLists):
  '''From a list of sequence chunks, build a list of those sequences fasta strings'''
  seqStrs = []
  for seqList in seqLists:
    fastaList = [f'>seq{i+1}\n{seq}\n' for i, seq in enumerate(seqList)]
    seqStrs.append(''.join(fastaList).strip())
  return seqStrs

def getFastaStrs(seqDic, maxChunk=1):
  '''Build a list of fasta strings from a list of sequences in chunks of maxChunk size'''
  maxChunk = len(seqDic) if not maxChunk else maxChunk
  seqList = list(seqDic.values())
  seqLists = divide_chunks(seqList, maxChunk)
  return buildSeqFasta(seqLists)

def getFastaFiles(seqDic, evalSoft, maxChunk=1):
  '''Write a series of fasta files with maxChunk number of sequences from a set of sequences'''
  fastaStrs = getFastaStrs(seqDic, maxChunk)
  faFiles = []
  for i, fStr in enumerate(fastaStrs):
    faFiles.append(f'/tmp/{evalSoft}_input_{i}.fa')
    with open(faFiles[-1], 'w') as f:
      f.write(fStr)
  return faFiles

def setData(driver, paramDic):
  '''Sets the additional data parameters in the web of the software evaluation
  driver: selenium driver, with url set in the software web
  paramDic: dic, containing the parameter names (as dic keys) and corresponding values (dic values)
  '''
  from selenium.webdriver.common.by import By
  for dk, dv in paramDic.items():
    dataElements = driver.find_elements(By.NAME, dk)
    for dEl in dataElements:
      if dEl.get_attribute('value') == dv:
        dEl.click()
  return driver


def getDriver(browserData):
  from selenium import webdriver
  from selenium.webdriver.chrome.options import Options as ChromeOptions
  from selenium.webdriver.firefox.options import Options as FireOptions
  '''Return a selenium WebDriver object, depending on the selected browser
  - browserData: dic, contains the information about the browser to be used
    - name: str, the name of the browser to use (either "Chrome" for Google-Chrome or Firefox)
    - path: str, path for the browser executable in case of non default
  '''
  if not 'name' in browserData or browserData['name'] != 'Firefox':
    options = ChromeOptions()
    driverObj = webdriver.Chrome
    browserPath = '/usr/bin/google-chrome' if (not 'path' in browserData or not browserData['path'])\
      else browserData['path']
  else:
    options = FireOptions()
    driverObj = webdriver.Firefox
    browserPath = '/usr/bin/firefox' if (not 'path' in browserData or not browserData['path']) \
      else browserData['path']

  options._binary_location = browserPath
  options.add_argument('--headless')
  options.add_argument('--remote-debugging-pipe')
  driver = driverObj(options=options)
  return driver


def performRequest(seqKeys, driver, softData):
  from selenium.webdriver.common.by import By
  '''Performs a request in a evaluation software using selenium to emulate the browser.
  - seqData: dic, contains the keys and values of the web elements to write, including the
             sequence in the format expected by the web (fasta file, fasta string or sequence string)
  - driver: selenium driver to use for the request
  - softData: dic, containing all the characteristics and info for the specific sofware web. Among others (key: value):
    - url: str, evaluation software url
    - params: dict, additional data arguments to fill in the web form
    - submitCSS: str, css selector to identify the submit button (e.g: "input[name='Submit']")
  - seqKeys: dic, if not None, specifies the web html name key and value to write the sequence name. e.g: {seqName: seq1}
  '''
  driver.get(softData['url'])

  for xKeyName, xKeyVal in seqKeys.items():
    extraElem = driver.find_element(By.NAME, xKeyName)
    extraElem.send_keys(xKeyVal)

  driver = setData(driver, softData['params'])
  driver.find_elements(By.CSS_SELECTOR, softData['submitCSS'])[0].click()
  return driver


def getSeqData(seqDic, softData):
  '''Returns a list containing the chunks of sequences as expected from the web to use.
  It can be either: a list with one fasta file, a list with one fasta string or a list with sequences strings
  - seqDic: dic, sequences {seqId: seqString}
  - softData: dic, containing all the characteristics and info for the specific sofware web. Among others (key: value):
    - multi: whether the web admits multiple sequences at one time
    - seqFormat: whether to return a fasta file ("fastaFile") or the fasta string ("fastaString")
    - softName: software name for the fasta file to be named
  '''
  if softData['multi']:
    if softData['seqFormat'] == 'fastaFile':
      seqData = getFastaFiles(seqDic, softData['softName'])
    else:
      seqData = getFastaStrs(seqDic)
  else:
    seqData = list(seqDic.values())
  return seqData


def updateBatchDic(outDic, batchDic):
  '''Updates(appends) the lists inside the outDic values with the ones in the batchDic
  '''
  for key, values in batchDic.items():
    if key in outDic:
      outDic[key] += values
    else:
      outDic[key] = values
  return outDic


def seleniumRequest(seqDic, softData, browserData, parseFunction, seqNameKey=None):
  '''Perform a series of Selenium requests an operations to emulate the evaluation of a set of sequences by a software
  web server.
  - seqDic: dic, sequences {seqId: seqString}
  - softData: dic, contains the information necessary to build the software web request
  - browserData: dic, contains the information necessary to build the Selenium driver
  - parseFunction: func, parses the driver data once the request is performed and returns a dic {'Score' [sc1, ...]}
  - seqNameKey: str, if not None, include the sequence name as a web element value to write in this key
  '''
  # url, data, softName, seqFormat='fastaString', seqName='sequence', multi=True
  driver = getDriver(browserData)
  seqData = getSeqData(seqDic, softData)

  # Performing one request for each chunk of admitted data (just once if fasta admitted)
  outDic = {}
  for i, seq in enumerate(seqData):
    curSeqKeys = {softData['seqName']: seq}
    if seqNameKey:
      curSeqKeys.update({seqNameKey: f'seq{i + 1}'})

    driver = performRequest(curSeqKeys, driver, softData)
    # Parse the driver with the corresponding function for each software
    batchDic = parseFunction(driver)
    outDic = updateBatchDic(outDic, batchDic)
  return outDic


def innerSplit(text, preText, endText):
  results, splitted = [], text.split(preText)[1:]
  for text in splitted:
    results.append(text.strip().split(endText)[0].strip())
  return results

########## REQUESTS ##########

def makeRequest(url, action='post', data={}, headers={}):
  if action == 'post':
    response = requests.post(url, data=data, headers=headers)
  else:
    response = requests.get(url, data=data, headers=headers)

  if response.status_code == 200:
    pass
  else:
    print(f"There was an error in request to {url}: {response.status_code}")
  return response


########### SELENIUM CALLS ################

def callABCpredSelenium(seqDic, browserData={}, data={}):
  data = {"window": "16", "filter": 'on', 'Threshold': "0.51"} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html",
              'multi': False,
              'seqName': 'SEQ', 'params': data, 'submitCSS': "input[value='Submit sequence']"}

  outDic = seleniumRequest(seqDic, softData, browserData, parseABCpred, seqNameKey='SEQNAME')
  return outDic

def callABCpred(protsDic, data={}):
  'https://webs.iiitd.edu.in/raghava/abcpred/ABC_submission.html'
  oriUrl = "https://webs.iiitd.edu.in"
  data = {"window": "16", "filter": 'on', 'Threshold': "0.51"} if not data else data
  headers = {"Referer": os.path.join(oriUrl, "raghava/abcpred/ABC_submission.html")}

  outDic = {}
  for seqId, sequence in protsDic.items():
    data.update({"SEQ": sequence})

    submitUrl = os.path.join(oriUrl, "cgibin/abcpred/test1_main.pl")
    response = makeRequest(submitUrl, 'post', data, headers)
    outDic[seqId] = getABCpredScore(parseABCpredOutHTML(response))
  return outDic

def callLBtope(sequences, browserData={}, data={}):
  data = {"for": 'flx'} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/lbtope/protein.php",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Submit antigen for prediction']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseLBtope)
  return outDic


def callToxinPred(sequences, browserData={}, data={}):
  data = {'method': '8', 'eval': '10', 'thval': '0.0'} if not data else data
  softData = {'url': "https://webs.iiitd.edu.in/raghava/toxinpred/multi_submit.php",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Run Analysis!']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseToxinPred)
  return outDic


def callToxinPred2(sequences, browserData={}, data={}):
  data = {'terminus': '4', 'svm_th': '0.6'} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/toxinpred2/batch.html",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Submit']"}

  # todo: check when sequences >=19
  outDic = seleniumRequest(sequences, softData, browserData, parseToxinPred2)
  return outDic


def callIFNepitope(sequences, browserData={}, data={}):
  data = {"method": 'svm'} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/ifnepitope/predict.php",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'sequence', 'params': data, 'submitCSS': "input[value='Submit Peptides for Prediction']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseIFNepitope)
  return outDic


def callIL4pred(sequences, browserData={}, data={}):
  data = {"method": '3'} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/il4pred/predict.php",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Virtual Screening']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseToxinPred)
  return outDic


def callIL10pred(sequences, browserData={}, data={}):
  data = {"method": '1'} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/il10pred/predict3.php",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Run Analysis!']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseIL10pred)
  return outDic


def callAlgPred2(sequences, browserData={}, data={}):
  data = {"terminus": '4', 'svm_th': "0.3"} if not data else data

  softData = {'url': "https://webs.iiitd.edu.in/raghava/algpred2/batch.html",
              'multi': True, 'seqFormat': 'fastaString',
              'seqName': 'seq', 'params': data, 'submitCSS': "input[value='Submit']"}

  outDic = seleniumRequest(sequences, softData, browserData, parseAlgPred2)
  return outDic


############## PARSING ##############

def getABCpredScore(resDic):
  '''Returns the dictionary in the same format as the rest of epitope selectors
  :param resDic: {'Rank': [], 'Sequence': [], 'Position': [], 'Score': []}
  :return: {(position, epitopeString): score}
  '''
  outDic = {}
  for i, epitope in enumerate(resDic['Sequence']):
    outDic[(resDic['Position'][i], epitope)] = resDic['Score'][i]
  return outDic


def parseABCpredOutHTML(response):
  '''Parse the ABCpred web server response table
  :param response: response of post to ABCpred server
  :return: {'Rank': [], 'Sequence': [], 'Position': [], 'Score': []}
  '''
  outDic = {}
  table = response.text.split('<table')[2]
  i, ks = 0, []
  for j, td in enumerate(table.split('<TD WIDTH')[1:]):
    value = td.split('>')[1].split('<')[0]
    if i == 0:
      ks.append(value)
      outDic[value] = []
    else:
      outDic[ks[j % len(ks)]].append(value)

    if 'TR>' in td:
      i += 1
  return outDic


def filterBestEpitopes(resDic, minProb=78):
  # todo: use input minprob
  epDic = {}
  for protId in resDic:
    epDic[protId] = {}
    for i, ep in enumerate(resDic[protId]):
      prob = resDic[protId][ep]['Probability']
      if prob >= minProb:
        epDic[protId] = updateBatchDic(epDic[protId],
                                       {'Sequence': [ep], 'Position': [i + 1], 'Score': [resDic[protId][ep]['Score']]})
  return epDic


def parseABCpred(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[width='60% bgcolor=']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[width='60% bgcolor=']")
  headerText = data[0].text
  seqName = innerSplit(headerText, 'Sequence name', '\n')[0]

  data = driver.find_elements(By.CSS_SELECTOR, "table[width='75% bgcolor=']")
  resultWeb = data[0]

  resDic = {}
  tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
  thead, tbody = tbody.find_elements(By.TAG_NAME, 'tr')[0], tbody.find_elements(By.TAG_NAME, 'tr')[1:]
  for cell in thead.find_elements(By.TAG_NAME, 'td'):
    resDic[cell.text] = []

  labels = list(resDic.keys())
  for row in tbody:
    for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
      if i < len(labels):
        resDic[labels[i]].append(cell.text)

  resDic['Position'] = resDic['Start position']
  del resDic['Start position']
  return {seqName: resDic}

def parseLBtope(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.PARTIAL_LINK_TEXT, 'Download results as a text file')
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.PARTIAL_LINK_TEXT, 'Download results as a text file')
  data[0].click()

  resTxt = driver.find_element(By.XPATH, "/html/body").text

  resDic = {}
  for line in resTxt.split('\n'):
    sline = line.split()
    if len(sline) == 3:
      ep, sc, perc = sline
      if 'X' not in ep:
        resDic[protId][ep] = {'Score': float(sc), 'Probability': float(perc)}
    else:
      protId = sline[2].replace('>', '')
      resDic[protId] = {}
  epDic = filterBestEpitopes(resDic)
  return epDic

def parseToxinPred11(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.ID, "tableTwo")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.ID, "tableTwo")
  resultWeb = data[0]

  resDic = {}
  thead = resultWeb.find_elements(By.TAG_NAME, 'thead')[0]
  for cell in thead.find_elements(By.TAG_NAME, 'b'):
    resDic[cell.text] = []

  labels = list(resDic.keys())
  tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
  for row in tbody.find_elements(By.TAG_NAME, 'tr'):
    for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
      resDic[labels[i]].append(cell.text)
  outDic = renameScore(resDic)
  return outDic

def parseToxinPred(driver):
  from selenium.webdriver.common.by import By
  '''Also used to parse IL4pred output since they use same template'''
  def getPages(driver):
    pageNum = driver.find_elements(By.CSS_SELECTOR, "input[class='pagedisplay']")[0]
    cPage, lastPage = pageNum.get_property('value').split('/')
    return cPage, lastPage

  def parseTable(resDic, resultWeb):
    labels = list(resDic.keys())
    tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
    for row in tbody.find_elements(By.TAG_NAME, 'tr'):
      for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
        resDic[labels[i]].append(cell.text)
    return outDic

  data = driver.find_elements(By.ID, "tableTwo")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.ID, "tableTwo")
  resultWeb = data[0]

  outDic = {}
  thead = resultWeb.find_elements(By.TAG_NAME, 'thead')[0]
  for cell in thead.find_elements(By.TAG_NAME, 'th'):
    outDic[cell.text] = []

  outDic = parseTable(outDic, resultWeb)
  cPage, lastPage = getPages(driver)
  while cPage != lastPage:
    nextPage = driver.find_elements(By.CSS_SELECTOR, "img[class='next']")
    driver.execute_script("arguments[0].click();", nextPage[0])
    outDic = parseTable(outDic, resultWeb)
    cPage, lastPage = getPages(driver)
  return renameScore(outDic)

def parseToxinPred2(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[border='1']")
  while not data:
    time.sleep(5)
    data = driver.driver.find_elements(By.CSS_SELECTOR, "table[border='1']")
  resultWeb = data[0]

  resDic = {}
  thead = resultWeb.find_elements(By.TAG_NAME, 'thead')[0]
  for cell in thead.find_elements(By.TAG_NAME, 'b'):
    resDic[cell.text] = []

  labels = list(resDic.keys())
  tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
  for row in tbody.find_elements(By.TAG_NAME, 'tr'):
    for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
      resDic[labels[i]].append(cell.text)
  outDic = renameScore(resDic)
  return outDic

def parseIFNepitope(driver):
  from selenium.webdriver.common.by import By
  def parseTable(outDic, resultWeb, labels):
    for row in resultWeb.find_elements(By.TAG_NAME, 'tr'):
      for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
        outDic[labels[i]].append(cell.text)
    return outDic

  data = driver.find_elements(By.ID, "example")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.ID, "example")
  resultWeb = data[0]

  labels = ['N0', 'Name', 'Epitope', 'Method', 'Result', 'Score']
  resDic = {l: [] for l in labels}

  outDic = parseTable(resDic, resultWeb, labels)
  nextPage = driver.find_elements(By.CSS_SELECTOR, "a[class='paginate_enabled_next']")
  while nextPage:
    nextPage[0].click()
    outDic = parseTable(outDic, resultWeb, labels)
    nextPage = driver.find_elements(By.CSS_SELECTOR, "a[class='paginate_enabled_next']")

  return outDic

def parseIL10pred(driver):
  from selenium.webdriver.common.by import By
  def parseTable(resDic, resultWeb):
    labels = list(resDic.keys())
    tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
    for row in tbody.find_elements(By.TAG_NAME, 'tr'):
      for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
        resDic[labels[i]].append(cell.text)
    return resDic

  data = driver.find_elements(By.CSS_SELECTOR, "table[class='table table-hover']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[class='table table-hover']")
  resultWeb = data[0]

  resDic = {}
  thead = resultWeb.find_elements(By.TAG_NAME, 'thead')[0]
  for cell in thead.find_elements(By.TAG_NAME, 'th'):
    resDic[cell.text] = []

  outDic = parseTable(resDic, resultWeb)
  nextPage = driver.find_elements(By.CSS_SELECTOR, "li[class='page-next']")
  while nextPage:
    nextPageBut = nextPage[0].find_elements(By.TAG_NAME, 'a')[0]
    nextPageBut.click()
    outDic = parseTable(outDic, resultWeb)
    nextPage = driver.find_elements(By.CSS_SELECTOR, "li[class='page-next']")

  return renameScore(outDic)

def parseAlgPred2(driver):
  from selenium.webdriver.common.by import By
  data = driver.find_elements(By.CSS_SELECTOR, "table[border='1']")
  while not data:
    time.sleep(5)
    data = driver.find_elements(By.CSS_SELECTOR, "table[border='1']")
  resultWeb = data[0]

  resDic = {}
  thead = resultWeb.find_elements(By.TAG_NAME, 'thead')[0]
  for cell in thead.find_elements(By.TAG_NAME, 'th'):
    resDic[cell.text] = []

  labels = list(resDic.keys())
  tbody = resultWeb.find_elements(By.TAG_NAME, 'tbody')[0]
  for row in tbody.find_elements(By.TAG_NAME, 'tr'):
    for i, cell in enumerate(row.find_elements(By.TAG_NAME, 'td')):
      resDic[labels[i]].append(cell.text)
  outDic = renameScore(resDic)
  return outDic

def renameScore(outDic, scoreKey=''):
  '''Rename the score key in a dict with just "Score"'''
  scoreK = None
  for k in outDic:
    if (scoreKey and scoreKey == k) or (not scoreKey and 'score' in k.lower()):
      scoreK = k
  outDic['Score'] = outDic.pop(scoreK)
  return outDic

def mapEvalParamNames(sDic):
  wsDic = {}
  for sName, curSDic in sDic.items():
    wsDic[sName] = {}
    for paramName, paramValue in curSDic.items():
      if paramName in EVAL_PARAM_MAP:
        paramName = EVAL_PARAM_MAP[paramName]
      wsDic[sName][paramName] = paramValue
  return wsDic