using ColorTypes, FileIO, ImageTransformations

function extract_image_from_file(fileName::String)
	if !isfile(fileName)
		println("""The file \""""*fileName*"""\" is missing.""")
		exit()
	else
		image = RGBA.(load(fileName))
		(r, g, b, a) = (red.(image), green.(image), blue.(image), alpha.(image))
		image = cat(r, g, b, a; dims = 3)
		image = round.(Int64, 255*image)
		return image::Array{Int64, 3}
	end
end

function overlay(images::Vector{Array{Int64, 3}})
	subBuffer = images[1][:, :, 1:3]
	for i in 2:length(images)
		subBuffer[:, :, 1] = round.(Int64, subBuffer[:, :, 1] + (images[i][:, :, 1] - subBuffer[:, :, 1]) .* (images[i][:, :, 4] / 255))
		subBuffer[:, :, 2] = round.(Int64, subBuffer[:, :, 2] + (images[i][:, :, 2] - subBuffer[:, :, 2]) .* (images[i][:, :, 4] / 255))
		subBuffer[:, :, 3] = round.(Int64, subBuffer[:, :, 3] + (images[i][:, :, 3] - subBuffer[:, :, 3]) .* (images[i][:, :, 4] / 255))
	end
	return subBuffer::Array{Int64, 3}
end

function deinterlace(image::Array{Int64, 3})
	deinterlaced = zeros(Int64, size(image))
	for i in 1:2:size(image)[1]
		deinterlaced[i, :, :] = image[i, :, :]
		if i+1 <= size(image)[1]
			deinterlaced[i+1, :, :] = image[i, :, :]
		end
	end
	return deinterlaced::Array{Int64, 3}
end

fileNames = [
"sample_131";
#"sample_132";
#"sample_133";
"sample_134";
#"sample_135";
#"sample_136";
#"sample_137";
"sample_138";
#"sample_139";
#"sample_140";
#"sample_141";
"sample_142";
#"sample_143";
#"sample_144";
#"sample_145";
"sample_146";
#"sample_147";
#"sample_148";
#"sample_149";
"sample_150"]
#"sample_151";
#"sample_152"]

scalebar = extract_image_from_file("100um_scalebar.png") # 640x640 px

images = [extract_image_from_file("data/"*fileName*".bmp") for fileName in fileNames]

for i in eachindex(images)
	
	#scaledImg = round.(Int64, imresize(images[i], (797, 640))[79:718, :, :])
	scaledImg = round.(Int64, imresize(deinterlace(images[i]), (797, 640))[79:718, :, :])
	scaledImgScalebar = overlay([scaledImg, scalebar])[:, :, 1] ./ 255
	
	save("output/"*fileNames[i]*"processed.png", scaledImgScalebar)
	
end
