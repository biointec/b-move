k = 4
with open(f"{k}/searches.txt", "r") as file:
    # Read each line in the file
    for line in file:
        # Split the line into three arrays
        array1, array2, array3 = line.split()

        # Remove the curly brackets and split the first array into elements
        elements = array1.strip("{}").split(",")

        # Convert the elements to integers
        elements = [int(element) for element in elements]

        # Mirror the order of the elements
        elements = [k - e for e in elements]

        # Convert the mirrored elements back to strings
        elements = [str(element) for element in elements]

        # Join the mirrored elements with commas
        mirrored_array1 = "{" + ",".join(elements) + "}"

        print(mirrored_array1, array2, array3)
