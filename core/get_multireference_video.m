% Author   : Philipp Flotho
% Copyright 2022 by Philipp Flotho, All rights reserved.

function video = get_multireference_video(video_reader, reference_idx, vid_number)
    video_start = cat(1, [1], find(diff(reference_idx) ~= 0));

    n_frames = length(reference_idx);
    if vid_number == length(video_start)
        idx = video_start(vid_number):n_frames;
    else
        idx = video_start(vid_number):video_start(vid_number + 1)-1;
    end
    video = SUBSET_file_reader(video_reader, idx);
end