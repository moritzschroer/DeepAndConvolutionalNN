import matplotlib.pyplot as plt
import numpy as np 

# Konvertiert das 1D-Array in ein 2D-Array, das zum Plotten der Ziffer verwendet wird.
def to2DArray(vektor):
    matrix = np.zeros((28,28))
    for i in range(28):
        for j in range(28):
            matrix[i, j] = vektor[i*28 + j]
            
    return matrix
    
# Plottet die Ziffer.
def plotDigit(bilddaten, verifizierungsdaten, geratene_ziffer = [None for i in range(9)]):
    fig, axs = plt.subplots(3, 3, figsize=(10, 10))

    for i in range(3):
        for j in range(3):
            axs[i, j].imshow(to2DArray(bilddaten[i*3 + j]), cmap='gray', interpolation='none')
            axs[i, j].set_title("digit: "+str(int(verifizierungsdaten[3*i + j]))+", prediction: "+str(geratene_ziffer[3*i + j]))
    plt.tight_layout()
    plt.show()

# Vorhersage von 9 zufällig ausgewählten Zahlen aus dem Testdatensatz mit dem Modell 
# und Plotten der Ziffern.
def predictAndPlot(modell, testdaten, test_verifikationsdaten):
    ziffern_verfügbar = range(len(testdaten))
    ziffern_zum_plotten = np.random.choice(ziffern_verfügbar, size = 9, replace=False)
    
    bilddaten, verifizierungsdaten, geratene_ziffer = [], [], []
        
    for j in range(9):
        vorhersage = modell.predict(np.array([testdaten[ziffern_zum_plotten[j]]]))
        bilddaten.append(testdaten[ziffern_zum_plotten[j]])
        verifizierungsdaten.append(test_verifikationsdaten[ziffern_zum_plotten[j]])
        geratene_ziffer.append(np.argmax(vorhersage))
        
    plotDigit(bilddaten, verifizierungsdaten, geratene_ziffer)