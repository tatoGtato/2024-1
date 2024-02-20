soker = imread("soker.jpg");
eye = imread("eye.jpg");


rect_eye = [111 33 65 58];
rect_soker = [163 47 143 151];
sub_eye = imcrop(eye,rect_eye);
sub_soker = imcrop(soker,rect_soker);


c = normxcorr2(sub_eye(:,:,1),sub_soker(:,:,1));
figure
surf(c) 
shading flat

% offset found by correlation
[max_c,imax] = max(abs(c(:)));
[ypeak,xpeak] = ind2sub(size(c),imax(1));
corr_offset = [(xpeak-size(sub_eye,2)) 
               (ypeak-size(sub_eye,1))];

% relative offset of position of subimages
rect_offset = [(rect_soker(1)-rect_eye(1)) 
               (rect_soker(2)-rect_eye(2))];

% total offset
offset = corr_offset + rect_offset;
xoffset = offset(1);
yoffset = offset(2);

xbegin = round(xoffset + 1);
xend   = round(xoffset + size(eye,2));
ybegin = round(yoffset + 1);
yend   = round(yoffset + size(eye,1));

% extract region from peppers and compare to onion
extracted_eye = soker(ybegin:yend,xbegin:xend,:);
if isequal(eye,extracted_eye) 
   disp("onion.png was extracted from peppers.png")
end

recovered_eye = uint8(zeros(size(soker)));
recovered_eye(ybegin:yend,xbegin:xend,:) = eye;
imshow(recovered_eye)