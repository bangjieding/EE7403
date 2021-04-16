from matplotlib import pyplot as plt
import cv2
import numpy as np

img_path = './figure/2.jpg'
img = cv2.imread(img_path)
width, height = img.shape[:2][::-1]

img_resize = cv2.resize(img, (int(width*0.5), int(height*0.5)), interpolation=cv2.INTER_CUBIC)
cv2.imshow("img", img_resize)
print("img_reisze shape:{}".format(np.shape(img_resize)))
# cv2.imwrite('./figure/img4_reshape.jpg', img_resize)

# img_gray = cv2.cvtColor(img_resize, cv2.COLOR_RGB2GRAY)
# cv2.imshow("img_gray",img_gray)
# print("img_gray shape:{}".format(np.shape(img_gray)))
# cv2.imwrite('./figure/img4_gray.jpg', img_gray)

bins = np.arange(257)
 
# item = img_gray
# hist,bins = np.histogram(item,bins)
# width = 0.7*(bins[1]-bins[0])
# center = (bins[:-1]+bins[1:])/2
# plt.bar(center, hist, align = 'center', width = width)
# plt.show()


chans = cv2.split(img_resize)
colors = ("blue","green","red")
for (chan,color) in zip(chans, colors):
    hist,bins = np.histogram(chan,bins)
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    plt.bar(center, hist, align = 'center', width = width, color=color)
plt.show()
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
cv2.waitKey()