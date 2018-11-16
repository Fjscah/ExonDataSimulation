import json
import re
from sys import argv

import matplotlib.pyplot as plt
import numpy
import pygal
from matplotlib.figure import Figure

from setting import READ


def char2phred(char,plus=33):
    code=ord(char)-plus
    return code

if __name__ == '__main__':
    filename ='SRR4238252_2.fastq'
    if len(argv)==2:
        script, filename= argv
    frequencies={}
    for x in range(1,READ+1):
        frequencies['pos%d_frequencies'%x]=[0]*43
    with open(filename,'r') as f:
        i =0
        for line in f.readlines():
            i+=1
            if i%4==0:
                row_fastq=line.strip()
                x=1
                for char in row_fastq:
                    frequencies['pos%d_frequencies'%x][char2phred(char)]+=1
                    x=x+1
            if i%100000==0:
                print(i)
    with open("frequencies.json","w") as f:
        json.dump(frequencies, f)
    '''
    def dividesum(sum):
        def divide(x):
            return x/sum
        return divide
    def addline(row):
        def addvalue(x):
            return x+row
        return addvalue

    with open("frequencies.json","r") as f:
        frequencies=json.load(f)

    for x in range(1,101):
        sum=numpy.sum(frequencies['pos%d_frequencies'%x])
        frequencies['pos%d_frequencies'%x]=list(map(dividesum(sum),frequencies['pos%d_frequencies'%x]))
        #frequencies['pos%d_frequencies'%x]=list(map(addline(100-x),frequencies['pos%d_frequencies'%x]))
        #print(numpy.sum(frequencies['pos%d_frequencies'%x]))
    x_values=[(10**(-x/10))/(1+10**(-x/10)) for x in range(43)]
    y_bottom=[0]*43
    ax=[]
    fig=[]

    for n in range(1):
        fig.append(plt.figure(figsize=(10, 150)))
        
        for x in range(n*100+1,n*100+101):
            ax.append(fig[n].add_subplot(100,1,(x-1)%100+1))
            ax[x-1].plot(x_values[0:20],frequencies['pos%d_frequencies'%x][0:20])
            ax[x-1].plot(x_values[0:20],y_bottom[0:20], c='blue', alpha=0.1)
            ax[x-1].fill_between(x_values[0:20],frequencies['pos%d_frequencies'%x][0:20],
            y_bottom[0:20], facecolor='blue', alpha=0.5)

    fig[0].savefig('qphred-percent_100.png', dpi=100)
    plt.show()
    '''
    '''
    plt.plot(x_values,frequencies['pos%d_frequencies'%1])
    plt.plot(x_values,y_bottom, c='blue', alpha=0.1)
    plt.fill_between(x_values,frequencies['pos%d_frequencies'%1],
        y_bottom, facecolor='blue', alpha=0.5)
    #plt.title("qphred of nucleaibase", fontsize=14)
    plt.xlabel("phred", fontsize=5)
    plt.ylabel("percent", fontsize=5)
    # 设置刻度标记的大小
    plt.tick_params(axis='both',  labelsize=5)
    '''
    '''
    chart = pygal.Bar(style=my_style, x_label_rotation=45, show_legend=False)
    chart.title = 'Most-Starred Python Projects on GitHub'
    chart.x_labels = names
    chart.add('', stars)
    chart.render_to_file('python_repos.svg')
    '''
