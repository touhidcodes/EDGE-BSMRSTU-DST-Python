#  fibonacci function
def fibonacci(n):
    if n <= 0:
        return "Input must be a positive integer."
    elif n == 1:
        return 0  
    elif n == 2:
        return 1  
    else:
        return fibonacci(n - 1) + fibonacci(n - 2)

# get fibonacci series
def generate_fibonacci_series(terms):
    if terms <= 0:
        print("Please enter a positive integer for the number of terms.")
        return
    print("Fibonacci Series:")
    for i in range(1, terms + 1):
        print(fibonacci(i), end=" ")
    print()

#  input number
num_terms = int(input("Enter the number of terms: "))
generate_fibonacci_series(num_terms)
