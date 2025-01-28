def main():
        # Get all FASTA files in the current directory
    fasta_files = [file for file in os.listdir() if file.endswith(".fasta")]
    
    #     # Get all files in the current directory
    # all_files = os.listdir()

    # # Filter for FASTA files
    # fasta_files = []
    # for file in all_files:
    #     if file.endswith(".fasta"):
    #         fasta_files.append(file)

    # Check if there are any FASTA files
    if not fasta_files:
        print("No FASTA files found in the current directory. Please add a .fasta file and try again.")
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
        if 0 <= file_index < len(fasta_files):
            selected_file = fasta_files[file_index]
            print(f"\nYou selected: {selected_file}")
        else:
            print("\nInvalid choice. Please select a valid number from the list.")
            return
    else:
        print("\nInvalid input. Please enter a number corresponding to the file.")
        return
