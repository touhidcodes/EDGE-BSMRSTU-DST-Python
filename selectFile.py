import os

def select_fasta_file():
    # Get all FASTA files in the current directory
    fasta_files = [file for file in os.listdir() if file.endswith(".fasta")]

    # Check if no FASTA files are found
    if not fasta_files:
        print("No FASTA files found in the current directory.")
        return
     
    # Display available FASTA files
    print("\nAvailable FASTA files:")
    for index, file_name in enumerate(fasta_files, 1):
        print(f"{index}. {file_name}")
    
      # Prompt user to select a file by number
    user_choice = input("\nEnter the number of the FASTA file you want to use: ")

    # Validate user input
    if user_choice.isdigit():
        file_index = int(user_choice) - 1
        if 0 <= file_index and file_index < len(fasta_files):
            selected_file = fasta_files[file_index]
            print(f"\nYou selected: {selected_file}")
            return selected_file
        else:
            print("\nInvalid choice. Please select a valid number from the list.")
            return
    else:
        print("\nInvalid input. Please enter a number corresponding to the file.")
        return


# Example usage
selected_file = select_fasta_file()
print(selected_file)
