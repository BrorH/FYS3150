from subprocess import run

start = "python master.py beam"
next = [
    "5 clear compile",
    "50",
    "50 plot=vec",
    "10,20,30,40,50 plot=vec",
    "10,20,30,40,50 plot",
    "5:101:5 plot=vec,counts",
]

description = [
    "n = 5, emtpy datafile, compile solver",
    "n = 50",
    "n = 50, plot smallest eigenvector",
    "n = 10,20,30,40,50, plot smallest eigenvectors for all",
    "same ns, plot smallest eigenvectors, all eigenvectors for n=50, and counts",
    "n = 5,10,15,...,90,95,100, plot smallest eigenvectors and counts",
]


for i in range(len(next)):
    print("description:")
    print(description[i])
    print("whats written:")
    print(f"{start} {next[i]}")
    print("output:")
    run(f"{start} {next[i]}".split())
    print("\n" * 3)
