import numpy

array = numpy.array([[2, 4, 6], [5, 10, 15], [6, 12, 18], [7, 14, 21], [8, 16, 24]])
print("Printing 2D Array")
print(array)

print("Choose 3 sample rows from 2D array")
randomRows = numpy.random.randint(100, size=30)
print(randomRows)
#for i in randomRows:
#    print(array[i, :])