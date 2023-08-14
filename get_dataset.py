import os
import fastaparser
import gc

dir_list = os.listdir("E:/VYSOKA/BAKA/data")
#dir_list = dir_list[306:]

def solve(nums):
    nums.sort()
    nums.append(1e9)
    ans=[]
    l=nums[0]
    for i in range(1,len(nums)):
        if nums[i] != nums[i-1] + 1:
            ans.append([l, nums[i-1]])
            l=nums[i]
    return ans


def get_data(dir_list):
    Done = 0  
    for dir in dir_list:
        file_list = os.listdir("E:/VYSOKA/BAKA/data/" + dir)
        Done = Done + 1
        print(Done,  "/491")

        #isExist = os.path.exists("E:/VYSOKA/BAKA/data/" + dir + "/" + dir +"_Offgene300.txt")
        #if isExist : continue

        #if "Geosmithia_morbida" in dir : break

        newFileList = []
        for file in file_list:
                if "genomic.gff.gz" in file : continue
                if "genomic.gff" in file : newFileList.append(file)
                if "genomic.fna.gz" in file : continue
                if "genomic.fna" in file : newFileList.append(file)

        newFileList.sort(key= lambda x: x.split(".")[1])
      
        path_to_fna = 'data/' + dir +'/' +newFileList[0]
        path_to_gff = 'data/' + dir +'/' +newFileList[1]
        
        with open(path_to_gff) as f:
            lines_gff = f.readlines()


        seq_regions = []

        CDS_in_region=[]
        for line in lines_gff:
            if "sequence-region" in line:
                seq_regions.append(CDS_in_region)
                CDS_in_region = []
            if "CDS" in line:
                l = line.split('\t')[3:9]
                l[5] = l[5].split(';')[0]
                l[0] = int(l[0])
                l[1] = int(l[1])
                CDS_in_region.append(l)


            if "###" in line.split(" ")[0]:
                seq_regions.append(CDS_in_region)

        seq_regions.pop(0)
        AllCDS = seq_regions

        def sortfucn(e):
            return e[0]

        for i in range(0,len(seq_regions)-1):
            if len(seq_regions[i]) ==0 :
                continue
            seq_regions[i].sort(key=sortfucn)

        AllCDS = seq_regions

        import copy
        chDuplicates = copy.deepcopy(seq_regions)
        for region in chDuplicates :
            for cds in region:
                chDuplicates[chDuplicates.index(region)][region.index(cds)] = cds[0]


        IndexesTodelete = []
        for region in chDuplicates:
            indexes = []
            newregion = []
            for cds in region:
                if cds not in newregion: 
                    newregion.append(cds)
                    continue
                index = [i for i, x in enumerate(region) if x == cds]
                for ind in index:
                    if ind in indexes:continue
                    indexes.append(ind)
            IndexesTodelete.append(indexes)

        for region in IndexesTodelete:
            deletedNumbers = 0
            for index in region:
                try:
                    del seq_regions[IndexesTodelete.index(region)][index- deletedNumbers]
                except:
                    deletedNumbers = deletedNumbers + 1
                    continue

        promotors = []

        for region in seq_regions:
            promotors_in_region = []
            for CDS in region:
                if CDS[3] == '+' :
                    promotors_in_region.append([int(CDS[0])-300,int(CDS[0])-1,"+",'n',CDS[5]])
                    continue
                promotors_in_region.append([int(CDS[1])+1,int(CDS[1])+300,"-",'n',CDS[5]])
            promotors.append(promotors_in_region)

        for region in promotors:
            for promotor in region:
                if promotor[2] == '+' :
                    if region.index(promotor) == 0 : continue
                    secondIndexOfCDSbefore = int(seq_regions[promotors.index(region)][region.index(promotor)-1][1])
                    if promotor[0] <= secondIndexOfCDSbefore :
                        promotor[3] = int(promotor[1]) - secondIndexOfCDSbefore
                        if promotor[3] < 0 :
                            promotor[3] = 0
                    continue
                if region.index(promotor)+1 == len(seq_regions[promotors.index(region)]) : continue
                firtsIndexOfCDSafter = int(seq_regions[promotors.index(region)][region.index(promotor)+1][0])
                if promotor[1] >= firtsIndexOfCDSafter :
                    promotor[3] = firtsIndexOfCDSafter - promotor[0]
                    if promotor[3] < 0 : promotor[3] = 0

        
        Numbers = []
        for region in seq_regions:
            nbr_in_region = []
            for numbers in region :
                number1 = int(numbers[0])
                number2 = int(numbers[1]) 
                nbr_in_region.append(list(range(number1,number2)))
            Numbers.append(nbr_in_region)

        
        for region in promotors :
            nbr_in_region = []
            for numbers in region :
                number1 = int(numbers[0])
                number2 = int(numbers[1]) 
                Numbers[promotors.index(region)].append((list(range(number1,number2))))

        
        Allnumbers = []
        for region in Numbers :
            wholestring = []
            for nbr in region:
                wholestring.append(str(nbr))
            string =''.join(wholestring)
            res = string.strip('][').replace('][',',').split(',')
            Allnumbers.append(res)

        


        Allnumbers = [[ subelt for subelt in elt if subelt != '' ] for elt in Allnumbers]

        
        IndexesOfOffgene = []
        for region in Allnumbers:
            res = [int(i) for i in region]
            IndexesOfOffgene.append(res)

        region_ranges = []
        offgene = []

        
        for line in lines_gff:
            if "sequence-region" in line:
                number1 =int(line.split(' ')[2])
                number2= int(line.split(' ')[3].replace("\n",''))
                region_ranges.append([number1,number2])

        for region in IndexesOfOffgene:
            rangeofregion = region_ranges[IndexesOfOffgene.index(region)]
            Allnumbers = list(range(rangeofregion[0],rangeofregion[1]))
            offgeneregion = list(set(Allnumbers) - set(region))
            offgene.append(offgeneregion)

        # PRIDAL
        OffgeneIntervals = []

        for region in offgene:
            OffgeneIntervals.append(solve(region))




        def Change_base(i):
            switcher = {
                "A": "T",
                "T": "A",
                "C": "G",
                "G": "C",
                "a": "t",
                "t": "a",
                "c": "g",
                "g": "c",
                }   
            return switcher.get(i,i)



        def rotate_seq(dataToRotate: str):
            reversedstring = dataToRotate[::-1]
            new_string = []
            for base in reversedstring :
                new_base = Change_base(base)
                new_string.append(new_base)
            new_data = "".join(new_string)
            return new_data

        FnaCDS = []
        FnaOff =[]
        Fnaprom = []
        with open(path_to_fna, "r") as rf:
            parser = fastaparser.Reader(rf, parse_method="quick")
            region = 0 
            for data in parser:
                header = data.header
                sequence = data.sequence
                for CDS in AllCDS[region] :
                    string = sequence[int(CDS[0])-1:int(CDS[1])]
                    if CDS[3] == "-":
                        string = rotate_seq(string)
                    FnaCDS.append([string,CDS[5],CDS[3]])
                for promotor in promotors[region]:
                    string = sequence[int(promotor[0])-1:int(promotor[1])]
                    if promotor[2] == "-":
                        string = rotate_seq(string)
                    Fnaprom.append([string,promotor[4],promotor[2],promotor[3]])
                for off in OffgeneIntervals[region] :
                    string = []
                    if off[1] - off[0] < 300 : continue
                    interval = list(range(off[0],off[1]))
                    
                    for offbase in interval: 
                        try:
                            string.append(sequence[offbase-1])
                        except:
                            continue
                    string = ''.join(string)
                    FnaOff.append([string,header.split(' ')[0]])
                region = region + 1

            
        pathtodata = "E:/VYSOKA/BAKA/data/"
        """ with open(pathtodata + dir + "/" + dir +"_promotor300.txt","w") as f:
            for prom in Fnaprom:
                f.write("\n" + "> " + prom[1]+ '   ' + prom[2]+ '   ' + str(prom[3]))
                line = prom[0]
                n = 80
                newline = [line[i:i+n] for i in range(0, len(line), n)]
                for line in newline:
                    f.write("\n" + line) """

        """ with open(pathtodata + dir + "/" + dir +"_CDS.txt","w") as f:
            for CDS in FnaCDS:
                f.write("\n" + CDS[1]+ '   ' + CDS[2])
                line = CDS[0]
                n = 80
                newline = [line[i:i+n] for i in range(0, len(line), n)]
                for line in newline:
                    f.write("\n" + line) """
        
         # PRIDAL TEST DO JMENA
        with open(pathtodata + dir + "/" + dir +"_Offgene300.txt","w") as f:
            for offs in FnaOff:
                f.write("\n" + offs[1])
                line = offs[0]
                n = 80
                newline = [line[i:i+n] for i in range(0, len(line), n)]
                for line in newline:
                    f.write("\n" + line)
        
        del promotors
        del Numbers
        del Allnumbers
        del IndexesOfOffgene
        del region_ranges
        gc.collect()

get_data(dir_list)

""" for index, dir in enumerate(dir_list):
    print(index, dir)
     """