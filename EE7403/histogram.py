from matplotlib import pyplot as plt
import cv2
import numpy as np

def HE_to_RGB(image):
    bins = np.arange(257)
    chans = cv2.split(image) #(b, g, r)
    colors = ("blue","green","red")
    for (chan, color) in zip(chans, colors):
        hist,bins = np.histogram(chan, bins)
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.bar(center, hist, align = 'center', width = width, color = color)
    plt.show()

    bH = cv2.equalizeHist(chans[0])
    gH = cv2.equalizeHist(chans[1])
    rH = cv2.equalizeHist(chans[2])

    chans2 = (bH, gH, rH)
    for (chan, color) in zip(chans2, colors):
        hist,bins = np.histogram(chan,bins)
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.bar(center, hist, align = 'center', width = width, color = color)
    plt.show()

    result = cv2.merge((bH, gH, rH))
    cv2.imshow("HE_RGB", result)
    cv2.imwrite('./figure/HE_RGB.jpg', result)

    cv2.waitKey()

def HE_to_HSV(image):
    imgHSV = cv2.cvtColor(image, cv2.COLOR_BGR2HSV)

    channelsHSV = cv2.split(imgHSV)
    channelsHSV[2] = cv2.equalizeHist(channelsHSV[2])

    channels = cv2.merge(channelsHSV)
    result = cv2.cvtColor(channels, cv2.COLOR_HSV2BGR)
    cv2.imshow("HE_HSV", result)
    cv2.imwrite('./figure/HE_HSV.jpg', result)

    bins = np.arange(257)
    chans = cv2.split(result) #(b, g, r)
    colors = ("blue","green","red")
    for (chan, color) in zip(chans, colors):
        hist,bins = np.histogram(chan, bins)
        width = 0.7*(bins[1]-bins[0])
        center = (bins[:-1]+bins[1:])/2
        plt.bar(center, hist, align = 'center', width = width, color = color)
    plt.show()


    cv2.waitKey()
    # imgYUV = cv2.cvtColor(image, cv2.COLOR_BGR2YCrCb)
    # cv2.imshow("src", image)

    # channelsYUV = cv2.split(imgYUV)
    # channelsYUV[0] = cv2.equalizeHist(channelsYUV[0])

    # channels = cv2.merge(channelsYUV)
    # result = cv2.cvtColor(channels, cv2.COLOR_YCrCb2BGR)
    # cv2.imshow("dst", result)

    # cv2.waitKey()


if __name__ == '__main__':
    img_path = './figure/8.jpg'
    img = cv2.imread(img_path)
    width, height = img.shape[:2][::-1]

    bins = np.arange(257)

    img_resize = cv2.resize(img, (int(width*0.5), int(height*0.5)), interpolation=cv2.INTER_CUBIC)
    cv2.imshow("img", img_resize)
    cv2.imwrite('./figure/8_resize.jpg', img_resize)

    img_gray = cv2.cvtColor(img_resize, cv2.COLOR_RGB2GRAY)
    cv2.imshow("img_gray",img_gray)
    cv2.imwrite('./figure/img8_gray.jpg', img_gray)


    img_eq = cv2.equalizeHist(img_gray)
    cv2.imwrite('./figure/img8_eq.jpg', img_eq)
    cv2.imshow("img_eq",img_eq)

    item = img_gray
    hist,bins = np.histogram(item,bins)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align = 'center', width = width)
    plt.show()
    
    item = img_eq
    hist,bins = np.histogram(item,bins)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align = 'center', width = width)
    plt.show()


    # HE_to_RGB(img_resize)
    # HE_to_HSV(img_resize)




# print("img_reisze shape:{}".format(np.shape(img_resize)))
# cv2.imwrite('./figure/img4_reshape.jpg', img_resize)

# img_gray = cv2.cvtColor(img_resize, cv2.COLOR_RGB2GRAY)
# cv2.imshow("img_gray",img_gray)
# print("img_gray shape:{}".format(np.shape(img_gray)))
# cv2.imwrite('./figure/img6_gray.jpg', img_gray)

# img_eq = cv2.equalizeHist(img_gray)
# # cv2.imwrite('./figure/img6_eq.jpg', img_eq)
# cv2.imshow("img_eq",img_eq)


# bins = np.arange(257)
 
# # normal = width * height * 0.25
# normal = 1

# item = img_gray
# hist,bins = np.histogram(item,bins)
# width = 0.7*(bins[1]-bins[0])
# center = (bins[:-1]+bins[1:])/2
# plt.bar(center, hist/normal, align = 'center', width = width)
# plt.show()

# item = img_eq
# hist,bins = np.histogram(item,bins)
# width = 0.7*(bins[1]-bins[0])
# center = (bins[:-1]+bins[1:])/2
# plt.bar(center, hist/normal, align = 'center', width = width)
# plt.show()


# chans = cv2.split(img_resize)
# colors = ("blue","green","red")
# for (chan,color) in zip(chans, colors):
#     hist,bins = np.histogram(chan,bins)
#     width = 0.7*(bins[1]-bins[0])
#     center = (bins[:-1]+bins[1:])/2
#     plt.bar(center, hist, align = 'center', width = width, color=color)
# plt.show()
# hist = cv2.calcHist([img_gray], [0], None,[256],[0,256])
# cv2.imshow('hist', hist)
# x = range(0,257)

# plt.figure()
# plt.title("Gray-scale Histogram")
# plt.xlabel("Bins")
# plt.ylabel("# of Pixels")
# plt.bar(x, hist)
# plt.xlim([0,256])
# plt.show()
# cv2.waitKey()

