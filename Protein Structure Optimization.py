import regex as re
import os
class Protein_Structure_optimization:
    def __init__(self,filename):
        self.filename=filename
    def initial_energy(self):
        os.system("./calRW "+self.filename+".pdb >> energy.txt")
        y=open("energy.txt", 'r')
        for data in y:
            list = re.split('[\s]+', data)
            y.close()
            return(abs(float(list[3])))
        y.close()
    def count_ca(self):
        file = open(self.filename+".pdb", 'r')
        count=0
        for data in file:
            if(re.search(r'^ATOM\s+\d+\s+CA\s+', data)):
                count+=1
        return(count)
    def write_file(self,i,new_point):
        file = open(self.filename+".pdb", 'r')
        j=1
        (x,y,z)=new_point
        s=""
        for data in file:
            if(re.search(r'^ATOM\s+\d+\s+CA\s+', data)):  
                if(i==j):
                    list = re.split('[\s]+', data)
                    f1,f2,f3=list[6],list[7],list[8]
                    x1,x2,x3=data.find(f1),data.find(f2),data.find(f3)
                    s+=data[:x1]+str(x)+data[x1+len(f1):x2]+str(y)+data[x2+len(f2):x3]+str(z)+data[x3+len(f3):]
                else:
                    s+=data
                j+=1
            else:
                s+=data
        y=open(self.filename+".pdb", 'w')
        y.write(s)
        y.close()
        os.system("./calRW "+self.filename+".pdb >> energy.txt")
        y=open("energy.txt", 'r')
        for data in y:
            list = re.split('[\s]+', data)
            energy=abs(float(list[3]))
            break
        y.close()
        return(energy)
    def read_file(self,i,delta):
        file = open(self.filename+".pdb", 'r')
        j=1
        s=["","","","","","",""]
        f1,f2,f3="0","0","0"
        for data in file:
            if(re.search(r'^ATOM\s+\d+\s+CA\s+', data)):  
                if(i==j):
                    list = re.split('[\s]+', data)
                    f1=f1[1:]+list[6]
                    f2=f2[1:]+list[7]
                    f3=f3[1:]+list[8]
                    x1,x2,x3=data.find(f1),data.find(f2),data.find(f3)
                    a1,a2,a3,a4,a5,a6=str(round(float(list[6])+(round(delta, 4)),4)),str(round(float(list[6])-(round(delta, 4)),4)),str(round(float(list[7])+(round(delta, 4)),4)),str(round(float(list[7])-(round(delta, 4)),4)),str(round(float(list[8])+(round(delta, 4)),4)),str(round(float(list[8])-(round(delta, 4)),4))
                    s[1]=s[0]+data[:x1]+a1+data[x1+len(f1):]
                    s[2]=s[0]+data[:x1]+a2+data[x1+len(f1):]
                    s[3]=s[0]+data[:x2]+a3+data[x2+len(f2):]
                    s[4]=s[0]+data[:x2]+a4+data[x2+len(f2):]
                    s[5]=s[0]+data[:x3]+a5+data[x3+len(f3):]
                    s[6]=s[0]+data[:x3]+a6+data[x3+len(f3):]
                    s[0]=""
                else:
                    s[0]=s[0]+data
                j+=1
            else:
                s[0]=s[0]+data
        s[1],s[2],s[3],s[4],s[5],s[6]=s[1]+s[0],s[2]+s[0],s[3]+s[0],s[4]+s[0],s[5]+s[0],s[6]+s[0]
        energynew=[]
        energynew.append(float(f1))
        energynew.append(float(f2))
        energynew.append(float(f3))
        for i in range(1,7):
            y=open(str(i)+".pdb", 'w')
            y.write(s[i])
            y.close()
            os.system("./calRW "+str(i)+".pdb >> energy.txt")
            y=open("energy.txt", 'r')
            for data in y:
                list = re.split('[\s]+', data)
                energynew.append(abs(float(list[3])))
                break
            y.close()
            y=open("energy.txt", 'w')
            y.write("")
            y.close()
        return(energynew)
    def solve(self,l):
        x=l[0]-(l[3]-l[4])//2
        y=l[1]-(l[5]-l[6])//2
        z=l[2]-(l[7]-l[8])//2
        return((x,y,z))
def main():
    delta=0.01
    epsilon=1
    filename="0"
    number_of_iteration=5
    obj=Protein_Structure_optimization(filename)
    energy=obj.initial_energy()
    print("intial",energy)
    count=obj.count_ca()
    print("count",count)
    new_energy=0
    for i in range(number_of_iteration):
        for j in range(1,count+1):
            l=obj.read_file(j,delta)
            new_energy=obj.write_file(j,obj.solve(l))
            print(new_energy)
        if(energy-new_energy<=epsilon):
            if(energy>new_energy):
                print("Final Energy = ",new_energy)
                return
            print("Final Energy = ",energy)
            return
        energy=new_energy
if __name__=="__main__":
    main()



    
