import numpy as np
import random
from jupyterNotebookProgressbar import log_progress


#this method reads a number of train_size handwritten digits and a number of
#test_size handwritten digits from the mnist dataset to arrays and returns the arrays

def readData(train_size, test_size):
    
    # Die Anzahl der handgeschriebenen Ziffern im MNIST-Datensatz beträgt 70000.
    if test_size + train_size > 70000:
        print("Fehler: Anzahl überschreitet Datenmenge")
        exit() 


    # Erzeugt ein numpy-Array mit Zahlen zwischen 0 und 59999.
    lines_available = np.arange(60000)

    # Wählt train_size zufällige Zeilenindizes aus lines_available aus, ohne Zurücklegen.
    train_lines_to_read = np.random.choice(lines_available, size=train_size, replace=False)

    # Entfernt die Indizes in train_lines_to_read aus lines_available.
    lines_available = np.setdiff1d(lines_available, train_lines_to_read)

    # Wählt test_size zufällige Indizes aus lines_available aus, ohne Zurücklegen.    
    test_lines_to_read = np.random.choice(lines_available, size=test_size, replace=False)

    # Nun werden die Zeilden mit den zufällig gewählten Indizes/Zeilennummern 
    # in die Arrays eingelesen.
    file = open('HandwrittenDigitData/mnist_data.csv', 'r')
    lines = file.readlines()
    file.close()

    train_data = np.zeros((train_size, 784))
    verify_train_data = np.zeros(train_size)

    print("Lese "+str(train_size)+" Trainingsdaten ein:")

    for j in log_progress(range(train_size)): 
            
        line = lines[train_lines_to_read[j]].split(',')
        verify_train_data[j] = np.int8(line[0])
        for k in range(784):
            train_data[j][k] = np.float32(line[k+1])
        
    test_data = np.zeros((test_size, 784))
    verify_test_data = np.zeros(test_size)

    print("Lese "+str(test_size)+" Testdaten ein")

    for j in log_progress(range(test_size)):
        line = lines[test_lines_to_read[j]].split(',')
        verify_test_data[j] = np.int8(line[0])
        for k in range(784):
            test_data[j][k] = np.float32(line[k+1])
        
    # Rückgabe der Arrays, wobei die Bilddatenarrays durch 256 geteilt werden, 
    # damit die Graustufen im Bereich [0,1] liegen
    return train_data/256, verify_train_data, test_data/256, verify_test_data

