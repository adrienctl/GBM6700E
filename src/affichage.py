import matplotlib.pyplot as plt

def plot_3D_points(mat3D):
    # CE SERAIT COOL DE FAIRE UN REPERE ORTHONORME
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for line in mat3D:
        ax.scatter(line[0],line[1],line[2],c='b',marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()