function Image(img)
  img.src = img.src:gsub('%.pdf$', '.svg')
  return img
end

