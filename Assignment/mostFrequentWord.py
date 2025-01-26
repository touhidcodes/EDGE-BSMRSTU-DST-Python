import os 

def find_most_frequent_word(file_name):
    # check if the file exists
    if os.path.exists(file_name):
        with open(file_name, 'r') as file:
            text = file.read().lower() 
            words = text.split()  
        
        max_frequency = 0
        max_frequency_word = ""
        
        # find the most frequent word
        for i in words:  
            count = 0
            for j in words: 
                if i == j:
                    count += 1
            if count > max_frequency:  
                max_frequency = count
                max_frequency_word = i
        
        print(f"The most frequent word is '{max_frequency_word}' and it appears {max_frequency} times.")
    else:
        print(f"The file '{file_name}' was not found. Please check the file name and try again.")


# input the file name
file_name = input("Enter the file name (with extension, e.g., 'text.txt'): ")
find_most_frequent_word(file_name)
