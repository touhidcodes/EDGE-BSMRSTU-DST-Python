def starPattern(num):
    
    #  star pattern
    for i in range(1, num + 1):
        print(" " * (num - i), end="")
        print("*" * i)

#  input number
pattern_num = int(input("Enter the number: "))
starPattern(pattern_num)
