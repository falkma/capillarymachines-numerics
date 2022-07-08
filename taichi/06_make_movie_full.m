%script to produce movie from saved images of taichi simulation
%adapted from https://www.mathworks.com/matlabcentral/answers/271490-create-a-movie-from-images

% load the images
movie_len = 400;
 images    = cell(2*movie_len,1);
 for i=1:movie_len
     images{i} = imread(sprintf('taichi_forward_movie_refine/taichi_movie%d.png',i));
 end
 for i=1:movie_len
     images{movie_len+i} = imread(sprintf('taichi_reverse_movie/taichi_movie%d.png',i));
 end
 % create the video writer with 10 fps
 writerObj = VideoWriter('taichi_full_movie.avi');
 writerObj.FrameRate = 10;
   % open the video writer
   open(writerObj);
   % write the frames to the video
    for u=1:(2*movie_len)
       % convert the image to a frame
       frame = im2frame(images{u});
       writeVideo(writerObj, frame);
   end
   % close the writer object
   close(writerObj);
   %this is the filename of the movie
   implay('taichi_full_movie.avi');