def nextmove(x, y):
    from random import randint
    """ randomly choose a direction
    0 = up, 1 = down, 2 = left, 3 = right"""
    direction =  randint(0,3)

    if direction == 0:  # move up
        y += 1
    elif direction == 1:  # move down
        y -= 1
    elif direction == 2:  # move right
        x += 1
    elif direction == 3:  # move left
        x -= 1
    else:
        print("error: direction isn't 0-3")

    return x, y

def Collision(x,y,grid):
    '''Checks for collision of particle
    with boundary or other particles'''
    import numpy as np
    if grid[y,x+1]==1 or grid[y,x-1]==1 or grid[y+1,x]==1 or grid[y-1,x]==1:
        return(True)
    else:
        return(False)