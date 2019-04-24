#list test script
import numpy as np

def mutate(transcript):
    listNew = transcript.copy()
    n = np.random.randint(0, len(listNew))
    listNew[n] = 'm'
    return listNew

i = 0
if i == 0:
        list1 = list(range(10))
        list2 = mutate(list1)
while i < 10:
    if i % 2 == 1:
        list1 = list2
        list2 = mutate(list1)
    else:
        list1 = list1
    print(list1, list2)
    i += 1
    