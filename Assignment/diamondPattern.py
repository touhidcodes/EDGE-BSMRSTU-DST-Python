def createDiamondPattern (num):
    
    #  for symmetry odd number is required
    if num % 2 == 0:
        print("You entered an even number. Don't worry, we'll adjust it to create an odd-numbered pattern for you!")
        num = num+ 1
    
    # upper triangle
    for i in range(num, 0, -2):
        print(" " * ((num - i) // 2), end="")
        print("*" * i)
    
    #  lower triangle
    for j in range(1, num+1, 2):
        print(" " * ((num - j) // 2), end="")
        print("*" * j)
        
createDiamondPattern(6);