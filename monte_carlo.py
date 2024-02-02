import numpy  as np
import matplotlib.pyplot as plt

def f(x,y,z):
    if abs(x) > np.pi/2 or abs(y) > np.pi/2 or abs(z) > np.pi/2 :
        return 0
    return np.cos(x*y*z)

def get_random_normal():
    return np.random.normal(),np.random.normal(),np.random.normal()

def get_normal_density_1d(x):
    return np.exp(-x**2/2)/(2*np.pi)**0.5

def get_normal_density_3d(x,y,z):
    return get_normal_density_1d(x) * get_normal_density_1d(y) * get_normal_density_1d(z)

def get_monte_carlo(f,n):
    S = 0
    for _ in range(n):
        x,y,z = get_random_normal()
        S += f(x,y,z)/get_normal_density_3d(x,y,z)
    return S/n

def plot_3d_function(f,n):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    for _ in range(n):
        x,y,z = get_random_normal()
        ax.scatter(x,y,z,c=(f(x,y,z)+1)/2,cmap="viridis")
    
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


if __name__ == "__main__":
    print(get_monte_carlo(f= f,n = 10000))
    plot_3d_function(f=f,n=100)
