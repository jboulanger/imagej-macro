# Tiles Tools

This is a macro to split an images into tiles, align the tiles and reassemble a montage.

Tiles are stored as "frames" (time points) in an hyperstack.

Hyperstack are supported for the splitting, apply and merge modes. For aligning multi-channel and multi-planes images, duplicate the stack to extract or combine channels and slices (using a maxiumum intensity projection for example), align the resulting stack and then apply the estimated shifts to the original hyperstack. 

The interface is the following:

![image](https://user-images.githubusercontent.com/3415561/123950857-0fe93780-d99c-11eb-9a46-a28fea8e7856.png)


## Modes
- Split: split the images into tiles with a determined overlap, the tiles are saved as frames in a new hyperstack
- Distribute: distribute the tiles spatially and allow to reposition the ROI
- Update: update the ROI position (x and y) in the result table
- Align: align the tiles from an hyperstack using the overlapping edges, shifts are saved in a table as columns dx and dy
- Apply: use the shift from the table to correct shifts for each frame of the hyperstack
- Merge: merge the frames into a mosaic fusing the overlaping tiles


