function Math(el)
  -- Replace all \bm{...} with \mathbf{...} inside math
  el.text = el.text:gsub("\\bm%s*{", "\\mathbf{")
  return el
end

