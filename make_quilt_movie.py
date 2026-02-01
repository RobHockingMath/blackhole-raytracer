from PIL import Image

# Define parameters
cols = 8
rows = 6

H = 560  # Height of each sub-image
W = 420  # Width of each sub-image

dN = round(768 / 2)
dn = 3  # Step size for selecting images
dn = 1

#for n in range(450+48,750-48):
for n in range(25,200-25):

    # Create a blank final image (black background)
    final_image = Image.new("RGB", (W * cols, H * rows), (0, 0, 0))

    count = 0

    # Iterate over rows and columns to place images
    for i in range(rows):
        print(f"Processing row {i+1}/{rows}")
        for j in range(cols):
            num = n-dn*24 + dn * count  # Compute image index

            # Generate filename
            #filename = f"./Oppenheimer_final/{num}.bmp"
            filename = f"./impossible_frames/{num:04d}.png"

            try:
                # Open image and resize using PIL
                image_ij = Image.open(filename).resize((W, H), Image.LANCZOS)#.transpose(Image.FLIP_TOP_BOTTOM)
                
                # Paste the image into the correct position in the final image
                final_image.paste(image_ij, (W * j, H * (rows - i - 1)))

            except FileNotFoundError:
                print(f"Warning: File {filename} not found. Skipping.")

            count += 1

    # Save the final image
    final_image.save(f"impossibleQuilt/impossible_{n-24}_im_qs8x6a0.75.png")

    print("Image saved successfully.")