from PIL import Image

def merge_images_with_gap(image1_path, image2_path, output_path, gap_width=250):
    # Open the two images
    img1 = Image.open(image1_path)
    img2 = Image.open(image2_path)

    # Flip both images vertically (reflection about the horizontal axis)
    img1 = img1.transpose(Image.FLIP_TOP_BOTTOM)
    img2 = img2.transpose(Image.FLIP_TOP_BOTTOM)

    # Ensure both images have the same height
    if img1.height != img2.height:
        new_height = min(img1.height, img2.height)
        img1 = img1.resize((int(img1.width * new_height / img1.height), new_height))
        img2 = img2.resize((int(img2.width * new_height / img2.height), new_height))

    # Create a new blank image with the combined width + gap and the same height
    new_width = img1.width + img2.width + gap_width
    merged_image = Image.new('RGB', (new_width, img1.height), color=(0, 0, 0))  # Black background

    # Paste the two images with a gap in between
    merged_image.paste(img1, (0, 0))
    merged_image.paste(img2, (img1.width + gap_width, 0))

    # Save the merged image
    merged_image.save(output_path)

# Process images from 0 to 99
for d in range(100):
    print(d)
    image1 = f'Oppenheimer_loop10/othereye_{d}.bmp'
    image2 = f'Oppenheimer_loop10/{d}.bmp'
    output = f'stereo/{d}.jpg'
    
    try:
        merge_images_with_gap(image1, image2, output, gap_width=250)
        print(f"Successfully merged: {output}")
    except Exception as e:
        print(f"Error processing {image1} and {image2}: {e}")