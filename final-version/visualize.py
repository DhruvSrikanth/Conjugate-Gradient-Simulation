import numpy as np
import matplotlib.pyplot as plt
import ast
import os

def read_file(filename):
    '''
    Read the file and return the data.
    '''
    with open(filename, "r") as f:
        file_str = f.read()
        data = ast.literal_eval(file_str)
        n_iter = data[0][0]
        grid = data[1]
        grid = np.array(grid)
    return (n_iter, grid)


def generate_plots(filename):
    '''
    Plot and save the data.
    '''
    n_iter, data = read_file(filename)
    plt.clf()
    plt.imshow(data)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.colorbar(orientation='vertical')
    plt.title('Conjugate Gradient Plot')

    plt.savefig('.' + filename.split('.')[1] + '.png')
    
    plt.close()

def generate_movie():
    '''
    Plot and save the data.
    '''
    imgs = []
    step = 100
    NT = len(os.listdir('./output/'))*step
    for i in range(0, NT, step):
        try:
            data = read_file("./output/output_x_" + str(i) + ".txt")
            imgs.append(data)
        except FileNotFoundError:
            continue
    
    fig = plt.figure()
    viewer = fig.add_subplot(111)
    plt.ion() # Turns interactive mode on (probably unnecessary)
    fig.show() # Initially shows the figure
    for i in range(len(imgs)):
        print("Plotting Simulation Sample : " + str(i+1))
        viewer.clear() # Clears the previous image
        viewer.imshow(imgs[i][1]) # Loads the new image
        plt.pause(.1) # Delay in seconds
        fig.canvas.draw() # Draws the image to the screen



if __name__ == '__main__':
    gen_mov = False

    # generate_plots('./output/mype_0.txt')
    # generate_plots('./output/mype_1.txt')
    generate_plots('./output/b.txt')
    generate_plots('./output/x.txt')
    

    if gen_mov:
        generate_movie()