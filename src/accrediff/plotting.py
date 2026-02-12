# accrediff/plotting.py
#import cv2 # type: ignore
#from natsort import natsorted # type: ignore
import os
#**************************************************************************************************************************************
def images_to_video(image_folder, output_video, frame_rate):
    try:
        import cv2
        from natsort import natsorted
    except ImportError:
        raise ImportError(
            "OpenCV is required for images_to_video()."
            " Please install it using 'pip install opencv-python'."
        )
    images = [img for img in os.listdir(image_folder) if img.endswith((".png", ".jpg", ".jpeg"))]
    #images.sort()  # Ensure the images are in the correct order
    images = natsorted(images) 

    # Read the first image to get the width and height
    first_image_path = os.path.join(image_folder, images[0])
    first_image = cv2.imread(first_image_path)
    height, width, layers = first_image.shape

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use 'XVID' or 'mp4v' for .mp4
    video = cv2.VideoWriter(output_video, fourcc, frame_rate, (width, height))

    count = 0
    for image in images:
        image_path = os.path.join(image_folder, image)
        img = cv2.imread(image_path)
        if img is not None and img.shape[0] == height and img.shape[1] == width:
            video.write(img)
            count += 1
        else:
            print(f"警告: 读取图片失败或尺寸不一致 {image_path}")

    # Release the video writer
    video.release()
    print(f"视频生成完成，写入帧数: {count}")
#**************************************************************************************************************************************