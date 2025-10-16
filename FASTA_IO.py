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

class FASTA_IO:
    
    import os, io, copy 
    
    RevCompDict={'A' : 'T',
                 'T' : 'A',
                 'C' : 'G',
                 'G' : 'C',
                 'R' : 'Y',
                 'Y' : 'R',
                 'S' : 'S',
                 'W' : 'W',
                 'K' : 'M',
                 'M' : 'K',
                 'B' : 'V',
                 'V' : 'B',
                 'D' : 'H',
                 'H' : 'D',
                 '-' : '-',
                 'N' : 'N',
                 ' ' : ' ',
                 'a' : 't',
                 't' : 'a',
                 'c' : 'g',
                 'g' : 'c',
                 'r' : 'y',
                 'y' : 'r',
                 's' : 'w',
                 'w' : 's',
                 'k' : 'm',
                 'm' : 'k',
                 'b' : 'v',
                 'v' : 'b',
                 'd' : 'h',
                 'h' : 'd',
                 'n' : 'n'}

    def __init__(self,path_to_file): # initializes the object creates connection to the file and set intnial values
        """
        FASTA_IO()
        
        inputs:
            path_to_file : path to indexed fasta reference file.
    
        """
        try: #try to open the file
            self.path_to_fasta_file=path_to_file
            self.path_to_fai_file=path_to_file+'.fai'
            if self.os.path.exists(self.path_to_fasta_file) and  self.os.path.exists(self.path_to_fai_file):
                self.fasta_file_handle = self.io.open(self.path_to_fasta_file,'r+')
                self.fai_file_handle = self.io.open(self.path_to_fai_file,'r')
            else:
                raise Exception
        except: #if file can not be opened create a file
            raise Exception("File or index could not be opened or does not exist")
        
        try: 
            self.__read_in_fai__()
        except:
            raise Exception("Could not process index")
    def __del__(self): #When object is deleted closes the connection to the file
        self.fasta_file_handle.close() 
        self.fai_file_handle.close() 

    def __eq__(self,other): # test if two IO objects are the same by comparing the file name will return true or false
        return self.path_to_files == other # compare file names

    def __ne__(self,other): # test if two IO objects are the different by comparing the file name will return true or false
        return not self.path_to_files == other # compare file names
    
    def __read_in_fai__(self):
        
        #######################################################################
        # Load in and process data from fai file
        #######################################################################
        fai_data = [line.strip().split('\t') for line in self.fai_file_handle.readlines()]
        self.SequenceIDs = [line[0].split(" ")[0] for line in fai_data] 
        self.Length = [int(line[1]) for line in fai_data] 
        self.Start = [int(line[2]) for line in fai_data] 
        self.BasesPerLine = [int(line[3]) for line in fai_data] 
        self.BytesPerLine = [int(line[4]) for line in fai_data] 
        self.fai_file_handle.seek(0,0)
        
    def get_sequence_IDs(self):
        """
        get all availble sequence IDs in the file

        get_sequence_IDs()

        Returns
        -------
        LIST 
            all avaialble sequence IDs found in the file .

        """
        return self.SequenceIDs
    def read_in_section(self,SequenceID,StartBase,StopBase):
        """
        read_in_section(SequenceID,StartBase,StopBase)
        
        inputs:
            SequenceID : the chromosome name e.g. 1 or X
            StartBase : start location of sequence to read (bed format) 
            StopBase : end location of sequence to read (bed format)
        outputs:
            returns the section from the start to the stop base in requested chromosome
        """
        
        if StartBase>StopBase:
            RevComp=True
            temp=self.copy.deepcopy(StartBase)
            StartBase=self.copy.deepcopy(StopBase)
            StopBase=self.copy.deepcopy(temp)
        else:
            RevComp=False
        
        try: #get the correct file handle for the chromsome
            SequenceIndex=self.SequenceIDs.index(SequenceID)
            File_handle=self.fasta_file_handle
            SequenceStart=self.Start[SequenceIndex]
            Length=self.Length[SequenceIndex]
            BasesPerLine=self.BasesPerLine[SequenceIndex]
            BytesPerLine=self.BytesPerLine[SequenceIndex]
            ByteTail=BytesPerLine-BasesPerLine
        except: #if file handle can not be found 
            raise Exception("incorrect chromosome id given")
            
        if StartBase<0:
            raise Exception("start position is less than 0")
            
        if StopBase>Length:
            raise Exception("stop position exceeds chromosome length")
            
            
            
            
        try: # get data start to stop
            number_of_new_line_at_start=(StartBase)//BasesPerLine
            number_of_new_line_at_stop=(StopBase)//BasesPerLine
            additional_bytes=(number_of_new_line_at_stop-number_of_new_line_at_start)*ByteTail
            length_to_read=StopBase-StartBase
            File_handle.seek(StartBase+SequenceStart+(number_of_new_line_at_start*ByteTail),0)
            data_out=File_handle.read(length_to_read+additional_bytes).replace('\n','').replace('\r','').replace(' ','')
            
            if RevComp:
                data_out="".join([self.RevCompDict.get(value,"N") for value in data_out[::-1]])
            
        except: #if could not get data
            raise Exception("Requested data outside the file range")
            
        return data_out
    
    def overwrite_section(self,SequenceID,StartBase,sequence):
        """
        overwrite_section(SequenceID,StartBase,sequence)
        
        inputs:
            SequenceID : the chromosome name e.g. 1 or X
            StartBase : start location of sequence to read (bed format) 
            sequence : sequence to write into file 
        outputs:
            returns the section from the start to the stop base in requested chromosome
        """
          
        if len(sequence)>0:
        
            StopBase=StartBase+len(sequence)
     
            try: #get the correct file handle for the chromsome
                SequenceIndex=self.SequenceIDs.index(SequenceID)
                File_handle=self.fasta_file_handle
                SequenceStart=self.Start[SequenceIndex]
                Length=self.Length[SequenceIndex]
                BasesPerLine=self.BasesPerLine[SequenceIndex]
                BytesPerLine=self.BytesPerLine[SequenceIndex]
                ByteTail=BytesPerLine-BasesPerLine
            except: #if file handle can not be found 
                raise Exception("incorrect chromosome id given")
                
            if StartBase<1:
                raise Exception("start position is less than 1")
                
            if StopBase>Length:
                raise Exception("stop position exceeds chromosome length")
                
                
            try: # get data start to stop
                number_of_new_line_at_start=(StartBase)//BasesPerLine
                File_handle.seek(StartBase+SequenceStart+(number_of_new_line_at_start*ByteTail),0)
                
                ### write the data for the line that contains the start position
                NumberBaseFirstLine=((number_of_new_line_at_start+1)*BasesPerLine)-StartBase
                first_chuck=sequence[:NumberBaseFirstLine]
                File_handle.write(first_chuck)
                
                ### write the data for each subseqent line
                sequence=sequence[NumberBaseFirstLine:]
                for chunk in [sequence[i:i+BasesPerLine] for i in range(0, len(sequence), BasesPerLine)]:
                    number_of_new_line_at_start+=1
                    File_handle.seek(SequenceStart+(number_of_new_line_at_start*BytesPerLine),0)
                    File_handle.write(chunk)
                    File_handle.flush()
                
            except: #if could not get data
                raise Exception("Could not write to file")

    def soft_mask_region(self,SequenceID,StartBase,StopBase):
        """
        soft_mask_region(SequenceID,StartBase,sequence)
        
        inputs:
            SequenceID : the chromosome name e.g. 1 or X
            StartBase : start location of sequence to read (bed format) 
            StopBase : end location of sequence to read (bed format)
        """
          
        
        self.overwrite_section(SequenceID,min(StartBase,StopBase),self.read_in_section(self,SequenceID,min(StartBase,StopBase),max(StartBase,StopBase)).lower())
        
        
    def hard_mask_region(self,SequenceID,StartBase,StopBase):
        """
        hard_mask_region(SequenceID,StartBase,sequence)
        
        inputs:
            SequenceID : the chromosome name e.g. 1 or X
            StartBase : start location of sequence to read (bed format) 
            StopBase : end location of sequence to read (bed format)
        """
          
        
        self.overwrite_section(SequenceID,StartBase,'N'*abs(StopBase-StartBase))
            