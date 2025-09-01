local struct = require 'pandoc.structure'

local fix_number = function (num)
  return num and num:gsub('^0%.', '') or nil
end

function Pandoc (doc)
  -- Create numbered sections
  doc.blocks = struct.make_sections(doc.blocks, {number_sections=true})

  -- Return the document with modified numbering
  return doc:walk {
    Div = function (div)
      local header = div.content[1]
      if not header or header.t ~= 'Header' then
        -- not a section div; return it unchanged
        return div
      end

      -- set the numbers attributes
      header.attributes.number = fix_number(header.attributes.number)
      div.attributes.number = fix_number(div.attributes.number)

      return div
    end
  }
end
