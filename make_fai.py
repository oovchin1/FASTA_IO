# -*- coding: utf-8 -*-
"""
MIT License

Copyright (c) 2025 Oleg S. Ovchinnikov

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""


def make_fai(path_to_file):
    """
    make_fai(path_to_file)
    
    inputs:
        path_to_file : path to indexed fasta reference file.

    """
    import io, os, copy
    
    if not os.path.exists(path_to_file):
        raise Exception("Fasta file does not exist")

    try: #try to open the file
        path_to_fasta_file=path_to_file
        path_to_fai_file=path_to_file+'.fai'
        Data_out=[]
        with io.open(path_to_fasta_file,'rb') as fasta_file_handle:
            
            line = fasta_file_handle.readline().decode('ascii')
            if line[:1] != ">":
                raise Exception()
                
            sequenceID=line[1:].split(' ')[0].replace('\n','').replace('\r','')
            sequenceStart=fasta_file_handle.tell()
            line = fasta_file_handle.readline().decode('ascii')
            BytesPerLine=len(line)
            line=line.strip()
            BasesPerLine=len(line)
            SequenceLength=len(line)
            NextLine = fasta_file_handle.readline().decode('ascii')
            
            while True:  
                line=copy.deepcopy(NextLine)
                NextLine=fasta_file_handle.readline().decode('ascii')
                
                if not NextLine: # Check if the line is empty (EOF)
                    SequenceLength+=len(line.strip())
                    Data_out.append([sequenceID,SequenceLength,sequenceStart,BasesPerLine,BytesPerLine])
                    break
                
                if line[:1] == ">":
                    Data_out.append([sequenceID,SequenceLength,sequenceStart,BasesPerLine,BytesPerLine])
                    
                    sequenceID=line[1:].split(' ')[0].replace('\n','').replace('\r','')
                    sequenceStart=fasta_file_handle.tell()
                    line = fasta_file_handle.readline().decode('ascii')
                    BytesPerLine=len(line)
                    line=line.strip()
                    BasesPerLine=len(line)
                    SequenceLength=len(line)
                    
                else:
                    if len(line) != BytesPerLine and NextLine[:1] != ">": 
                        raise Exception
                    line=line.strip()
                    if BasesPerLine != len(line) and NextLine[:1] != ">":
                        raise Exception
                    SequenceLength+=len(line)
                            
                
        with io.open(path_to_fai_file,'w') as fai_file_handle:
            fai_file_handle.write('\n'.join(['\t'.join([str(item) for item in value]) for value in Data_out]))

    except: #if file can not be opened create a file
        raise Exception("could not process Fasta File")
        raise Exception("File or index could not be opened or does not exist")
    